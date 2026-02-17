# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# (GRM_for_linear_LC)=
# # Effects of the Peclet number on the elution of a rectangular pulse 

# %% [markdown]
# This case study is a reproduction of parts of the results published in:
# <br><br>
# **Qamar et al. (2014)**:  
# *Analytical solutions and moment analysis of general rate model for linear liquid chromatography*,  
# *Chemical Engineering Science*, 107, 192–205.  
# [https://doi.org/10.1016/j.ces.2013.12.019](https://doi.org/10.1016/j.ces.2013.12.019)
# <br><br>
# This work on the General Rate Model (GRM) for linear liquid chromatography presented in that work was later adopted as a model problem to determine convergence benchmarks by:
# <br><br>
# **Leweke et al. (2018)**:  
# *Chromatography Analysis and Design Toolkit (CADET)*,  
# *Computers & Chemical Engineering*, 113, 274–294.  
# [https://doi.org/10.1016/j.compchemeng.2018.02.025](https://doi.org/10.1016/j.compchemeng.2018.02.025)
# <br><br>
#

# %% [markdown]
# In the given example, a tracer is introduced into the column as a **rectangular pulse**. The binding behavior follows the Linear binding model. The unit operation model for this process is the **General Rate Model**. In their study, Qamar et al. examined the influence of dispersion on the elution behavior by analysing the effects of different Peclet numbers (Fig.5). The **Peclet number `Pe`** is a value discribing the **ratio of convection to dispersion** for particle flow in a column. 
#
# Here, the elution is simulated for the numerical values from Table 2 (Leweke et al.). The interstitial velocity used here is u = 0.3 cm / min. Following this, the effects of different Peclet numbers are examined by simulating elutions for six different values.  
# The `flow_rate` can be calculated as the product of the interstitial cross section area and the interstitial velocity u. 

# %%
import numpy as np

from CADETProcess.processModel import ComponentSystem
from CADETProcess.processModel import Linear
from CADETProcess.processModel import Inlet, GeneralRateModel, Outlet
from CADETProcess.processModel import FlowSheet
from CADETProcess.processModel import Process

# Component System
component_system = ComponentSystem()
component_system.add_component('A')

# Binding Model
binding_model = Linear(component_system, name='Linear')
binding_model.is_kinetic = False
binding_model.adsorption_rate = [2.5]  # k_a [m_MP³ / (m_SP³ * s)]
binding_model.desorption_rate = [1.0]  # k_d [1 / s]


# Unit Operations
column = GeneralRateModel(component_system, name='column')
column.binding_model = binding_model
column.length = 0.017  # L [m]
column.cross_section_area = 1e-3 # [m²]
column.bed_porosity = 0.4  # ε_c [-]
column.particle_radius = 4.0e-5  # r_p [m]
column.particle_porosity = 0.333  # ε_p [-] 
column.axial_dispersion = 3.33e-9  # D_c [m² / s]
column.film_diffusion = column.n_comp * [1.67e-6]  # k_f [m / s]
column.pore_diffusion = column.n_comp * [3.003e-6]  # D_p [m² / s]
column.surface_diffusion = column.n_bound_states * [0.0]  # D_s [m² / s]
column.c = [0.0]  # [mM] 
column.q = [0.0]  # [mM]  

inlet = Inlet(component_system, name='inlet')
inlet.flow_rate = column.cross_section_area_interstitial * (0.3 * (1e-2 / 60))  # m² * [m / s] 

outlet = Outlet(component_system, name='outlet')


# Flow Sheet
flow_sheet = FlowSheet(component_system)

flow_sheet.add_unit(inlet)
flow_sheet.add_unit(column)
flow_sheet.add_unit(outlet, product_outlet=True)

flow_sheet.add_connection(inlet, column)
flow_sheet.add_connection(column, outlet)


# %%
# Process
process = Process(flow_sheet, 'pulse')
process.cycle_time = 100 * 60  # [s]
pulse_duration = 20.0 * 60  # [s]

c_pulse = np.array([1.0])  # injection concentration = 1.0 [mol / m³]
c_initial = np.array([0.0])
  
process.add_event('pulse_start', 'flow_sheet.inlet.c', c_pulse)
process.add_event('pulse_stop', 'flow_sheet.inlet.c',  c_initial, pulse_duration)


# %%
print(__name__)
if __name__ == '__main__':
    from CADETProcess.simulator import Cadet
    process_simulator = Cadet()

    simulation_results = process_simulator.simulate(process)
    simulation_results.solution.column.outlet.plot()

# %% [markdown]
# The Peclet number is given by the product of the column length and the interstiatial velocity devided by the dispersion. The `peclet_number` of the column used by Qamar et al. and shown in the plot above is around 255. 

# %%
column.axial_dispersion

# %%
peclet_number = (column.length * (0.3 * (1e-2 / 60))) / column.axial_dispersion[0]
peclet_number

# %% [markdown]
# The effects of different Peclet numbers on the elution curve can be examined by varying the dispersion rates within the column. The results show an approaching of the elution profile to a **rectangular output** as the **dispersion is reduced** and the **Peclet number increases**. The greater influence of **advective transport** in high Peclet numbers results in an elution that is more similar to the concentration profile at the inlet. Throughout the column the velocities of the particles within the mobile phase differ less with increasing Peclet numbers. The elution approaches the behavior of an **ideal plug flow reactor**. <br>
# Smaller Peclet numbers indicate a greater influence of dispersion. This results in a more gradual elution with larger retention times as the differences between particle velocities in the mobile phase increase. The resulting elution profile approaches the behavior of a **CSTR**.
#

# %%
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from CADETProcess.simulator import Cadet
    process_simulator = Cadet()
       
    fig, axes = plt.subplots(3,2, figsize = [10, 12]) 
    axes = axes.flatten()
 
    peclet_numbers = [1, 5, 25, 50, 100, 255]
    for i, peclet_number in enumerate(peclet_numbers):
        column.axial_dispersion = (column.length * (0.3 * (1e-2 / 60))) / peclet_number
        simulation_results = process_simulator.simulate(process)
        simulation_results.solution.column.outlet.plot(ax = axes[i])
        axes[i].set_title(f"Peclet number: {peclet_number}")

    plt.tight_layout()

# %%
