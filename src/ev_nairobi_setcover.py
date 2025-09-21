#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, math, warnings, numpy as np, pandas as pd, geopandas as gpd, networkx as nx, pyomo.environ as pyo, folium, osmnx as ox
from shapely.geometry import Point
warnings.filterwarnings("ignore", category=UserWarning)
ox.settings.log_console = True
ox.settings.use_cache = True
PLACE_NAME = "Nairobi, Kenya"
TIME_THRESHOLDS_MIN = [10, 15, 20]
HEX_SIZE_KM = 1.2
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

def km_to_deg(km): return km / 111.32
def build_hex_grid(boundary_gdf, hex_km):
    poly = boundary_gdf.geometry.unary_union
    minx, miny, maxx, maxy = poly.bounds
    dx = km_to_deg(hex_km); dy = dx * 3**0.5 / 2
    xs = np.arange(minx - dx, maxx + dx, dx); ys = np.arange(miny - dy, maxy + dy, dy)
    pts = []
    for j, y in enumerate(ys):
        off = 0 if j % 2 == 0 else dx/2
        for x in xs: pts.append(Point(x+off, y))
    g = gpd.GeoDataFrame(geometry=pts, crs="EPSG:4326")
    g = g[g.within(poly)].reset_index(drop=True)
    g["id"] = [f"d{idx:04d}" for idx in range(len(g))]; g["weight"] = 1.0
    return g[["id","weight","geometry"]]

def nearest_nodes(G, gdf): return ox.distance.nearest_nodes(G, gdf.geometry.x, gdf.geometry.y)

def compute_reach_sets(G, cand_nodes, demand_nodes, thresholds_min):
    node_to_dem = {}
    for j, dn in enumerate(demand_nodes): node_to_dem.setdefault(int(dn), []).append(j)
    cover = {}
    for tmin in thresholds_min:
        cutoff = tmin*60; cover[tmin] = {j: [] for j in range(len(demand_nodes))}
        for i, cn in enumerate(cand_nodes):
            lengths = nx.single_source_dijkstra_path_length(G, cn, weight="travel_time", cutoff=cutoff)
            for node_id in lengths:
                if node_id in node_to_dem:
                    for j in node_to_dem[node_id]: cover[tmin][j].append(i)
    return cover

def solve_set_cover(num_sites, num_demands, cover_dict, solvers=("highs","cbc","glpk")):
    m = pyo.ConcreteModel()
    m.I = pyo.RangeSet(0, num_sites-1); m.J = pyo.RangeSet(0, num_demands-1)
    m.y = pyo.Var(m.I, domain=pyo.Binary)
    def cover_rule(m, j):
        L = cover_dict[j]
        return (sum(m.y[i] for i in L) >= 1) if L else pyo.Constraint.Skip
    m.c = pyo.Constraint(m.J, rule=cover_rule)
    m.o = pyo.Objective(expr=sum(m.y[i] for i in m.I), sense=pyo.minimize)
    used = None; solver = None
    for s in solvers:
        try:
            solver = pyo.SolverFactory(s)
            if solver.available(False): used = s; break
        except: pass
    if not used: raise RuntimeError("No MILP solver found. Install highspy or coincbc or glpk.")
    solver.solve(m, tee=False)
    picks = [i for i in range(num_sites) if pyo.value(m.y[i])>0.5]
    uncovered = [j for j in range(num_demands) if not cover_dict[j]]
    return picks, uncovered, used

def make_map(boundary, demand, candidates, scenario_results, html_path):
    center = [boundary.geometry.centroid.y.values[0], boundary.geometry.centroid.x.values[0]]
    m = folium.Map(location=center, zoom_start=11, control_scale=True)
    folium.GeoJson(boundary.__geo_interface__, name="Nairobi Boundary",
                   style_function=lambda x: {"fillOpacity": 0, "color": "#333", "weight": 2}).add_to(m)
    for tmin, pack in scenario_results.items():
        picks = pack["selected_indices"]; covered_js = set(pack["covered_demands"])
        sel = folium.FeatureGroup(name=f"Selected Sites ({tmin} min)", show=(tmin==10))
        cov = folium.FeatureGroup(name=f"Covered Demand ({tmin} min)", show=(tmin==10))
        unc = folium.FeatureGroup(name=f"Uncovered Demand ({tmin} min)", show=False)
        for i in picks:
            rr = candidates.iloc[i]
            folium.Marker([rr.geometry.y, rr.geometry.x], icon=folium.Icon(icon="bolt", prefix="fa"),
                          popup=str(rr.get("name", rr["id"]))).add_to(sel)
        for j, rr in demand.reset_index(drop=True).iterrows():
            layer = cov if j in covered_js else unc
            folium.CircleMarker([rr.geometry.y, rr.geometry.x], radius=2, fill=True, fill_opacity=0.7).add_to(layer)
        sel.add_to(m); cov.add_to(m); unc.add_to(m)
    folium.LayerControl(collapsed=False).add_to(m)
    m.save(html_path)

def main():
    boundary = ox.geocode_to_gdf(PLACE_NAME).to_crs("EPSG:4326")
    G = ox.graph_from_polygon(boundary.geometry.unary_union, network_type="drive", retain_all=True, simplify=True)
    G = ox.add_edge_speeds(G); G = ox.add_edge_travel_times(G)
    # Candidates: malls & supermarkets
    tags = [{"shop":"mall"},{"shop":"supermarket"},{"amenity":"supermarket"}]
    frames = []
    for tag in tags:
        g = ox.features_from_polygon(boundary.geometry.unary_union, tags=tag)
        if g is not None and len(g)>0: frames.append(g)
    candidates = (pd.concat(frames, ignore_index=True).to_crs("EPSG:4326"))
    candidates["geometry"] = candidates.geometry.centroid
    candidates["id"] = [f"s{idx:04d}" for idx in range(len(candidates))]
    demand = build_hex_grid(boundary, HEX_SIZE_KM)
    # Save raw
    candidates.to_file(os.path.join(RESULTS_DIR, "candidates.geojson"), driver="GeoJSON")
    demand.to_file(os.path.join(RESULTS_DIR, "demand.geojson"), driver="GeoJSON")
    # Nodes & reachability
    cand_nodes = nearest_nodes(G, candidates); dem_nodes = nearest_nodes(G, demand)
    cover_maps = compute_reach_sets(G, cand_nodes, dem_nodes, TIME_THRESHOLDS_MIN)
    scenario_results = {}; rows = []
    for tmin in TIME_THRESHOLDS_MIN:
        picks, uncovered, used = solve_set_cover(len(candidates), len(demand), cover_maps[tmin])
        sel_df = candidates.iloc[picks].copy(); sel_df["scenario_min"] = tmin
        sel_df.drop(columns="geometry").to_csv(os.path.join(RESULTS_DIR, f"selected_sites_{tmin}min.csv"), index=False)
        covered_demands = sorted(set([j for j in range(len(demand)) if cover_maps[tmin][j]]))
        scenario_results[tmin] = {"selected_indices": picks, "uncovered_demands": uncovered,
                                  "covered_demands": covered_demands, "solver": used}
        rows.append({"scenario_min": tmin, "num_candidates_total": len(candidates),
                     "num_demand_total": len(demand), "selected_sites": len(picks),
                     "coverable_demands": len(covered_demands), "uncovered_demands": len(uncovered), "solver": used})
    pd.DataFrame(rows).to_csv(os.path.join(RESULTS_DIR, "coverage_summary.csv"), index=False)
    make_map(boundary, demand, candidates, scenario_results, os.path.join(RESULTS_DIR, "nairobi_charging_map.html"))
    print("Done. See results/ for outputs.")
if __name__ == "__main__":
    main()
