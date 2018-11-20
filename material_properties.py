# --- DEFINE MATERIAL PROPERTIES ----

PMMA = {"k": 0.0002,      # kW/mK [Vermesi_PhD]
        "rho": 1190,      # kg/m3 [Vermesi_PhD]
        "c": 1.606,       # kJ/kgK [Vermesi_PhD]
        "emissivity": 1,  # (-)
        "absorptivity": 1,  # (-)
        }

PA6 = {"k": 1,
       "rho": 1,
       "c": 1,
       "emissivity": 1,
       "absorptivity": 1,
       }

TIMBER = {"k": 1,
          "rho": 1,
          "c": 1,
          "emissivity": 1,
          "absorptivity": 1,
          }

ALUMINIUM = {"k": 0.167,       # kW/mK [Vermesi_PhD]
             "rho": 2700,      # kg/m3 [Vermesi_PhD]
             "c": 0.896,       # kJ/kgK [Vermesi_PhD]
             "m": 0.4,         # kg
             "emissivity": 1,  # (-)
             "absorptivity": 1,  # (-)
             }
