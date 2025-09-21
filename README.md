# Nairobi EV Charging — Set Cover (Malls & Supermarkets)

This project finds the **minimum number of EV charging stations** to cover demand across **Nairobi** within a **road travel-time threshold** (10/15/20 minutes). Candidate locations are **malls & large supermarkets** (parking + dwell time).

We use OpenStreetMap via **OSMnx** to pull the city boundary, road network, and POIs, build a demand hex-grid, compute reachability by travel-time, and solve a **Set Cover** model with **Pyomo** using a free solver (HiGHS/CBC/GLPK). Finally, we render an interactive **Folium** map of chosen sites and coverage.

## Repo Structure
