from tkinter import *
from tkinter import filedialog, messagebox
import csv

####################### READ TABLE ##############################

# Import a CSV file containing a 2-column table (e.g., temperature and property)
def import_csv():
    file = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
    if not file:
        return None
    try:
        with open(file, newline='') as f:
            reader = csv.reader(f, delimiter=';')
            next(reader)  # skip title
            data = [row for row in reader]
            return data
    except Exception as e:
        messagebox.showerror("Error", f"Error while importing:\n{e}")
        return None

# Read material table
def read_therm_table(data_therm, porosity, solidus, liquidus, ambient_temp):
    global data_density, data_density_powder, data_spe_heat, data_cond, data_cond_powder, ambient_density
    data_density = []
    data_density_powder = []
    data_cond_powder = []
    data_density_solid = []
    data_cond = []
    data_spe_heat = []
    ambient_density = None

    for row in data_therm:
        row = [value.replace(',', '.') for value in row[:4]]
        temp, density, conductivity, spe_heat = map(float, row)

        data_density.append((density, temp))
        if temp >= ambient_temp and ambient_density == None: # density for mechanical analysis
            ambient_density = density

        if temp == solidus:  # To store densitypow(T_sol) value
            density_sol = density
            conductivity_sol = conductivity
            spe_heat_sol = spe_heat

        # Compute mushy zone density (Scheil model)
        if temp == liquidus:
            delta_T = liquidus - solidus
            for i in range(1, 9):
                temperature = solidus + delta_T / 9 * i
                liq = (temperature ** 2 - 2 * solidus * temperature + solidus ** 2) / (delta_T ** 2)  # liquid proportion
                data_density_powder.append((round(density * liq + density_sol*(1 - porosity)*(1-liq), 2), round(temperature, 2), 0.0))  # density = (densitysol(T_liq)-densitypow(T_sol))*liq+densitypow(T_sol)
                data_density_solid.append((round((density - density_sol) * liq + density_sol, 2), round(temperature, 2), 1.0))
                data_cond_powder.append((round(conductivity* liq + conductivity_sol*(1 - porosity)*(1-liq), 2), round(temperature, 2))) 
                data_cond.append((round((conductivity - conductivity_sol) * liq + conductivity_sol, 2), round(temperature, 2)))
                data_spe_heat.append((round((spe_heat - spe_heat_sol) * liq + spe_heat_sol, 2), round(temperature, 2)))
            data_density_powder.append((density, temp, 0.0))
            data_cond_powder.append((conductivity, temp))

        if temp <= solidus:
            data_density_powder.append((round(density * (1 - porosity),2), temp, 0.0))
            data_cond_powder.append((round(conductivity * (1 - porosity),2), temp))
        
        data_density_solid.append((density, temp, 1.0))
        data_spe_heat.append((spe_heat, temp))
        data_cond.append((conductivity, temp))

    data_density_powder += data_density_solid  # Concatenate the 2 tables, 0.0 first

# Read mechanical properties table
def read_mecha_table(data_mecha):
    global data_expansion, data_elastic
    data_expansion = []
    data_elastic = []

    for row in data_mecha:
        row = [value.replace(',', '.') for value in row[:4]]
        temp, expan, young, poisson = map(float, row)
        data_expansion.append((expan, temp))
        data_elastic.append((young, poisson, temp))

# Read stress-strain table
def read_pla_table(data_pla):
    global data_plastic
    data_plastic = []
    strain_values = [float(s.replace(',', '.')) * 0.01 for s in data_pla[0][1:]]

    for row in data_pla[1:]:
        row = [cell.replace(',', '.') if cell != '' else '' for cell in row]
        temp = float(row[0])
        stress_values = [float(val) if val != '' else None for val in row[1:]]

        for strain, stress in zip(strain_values, stress_values):
            if stress is not None:
                data_plastic.append((stress, strain, temp))

    data_plastic.sort(key=lambda x: x[2])


###################### WRITE FILES ###############################

# Write constant.py for Python script
def write_constant_py(job, material, latent_heat, solidus, liquidus, entries_dflux):
    # Job
    lines = [f"JOB = '{job.get()}' # Name of the thermal job"]
    
    # Write scalar parameters
    for key in entries_dflux.keys():
        lines.append(f"{key} = {entries_dflux[key].get()}")

    # Computation
    lines.append("LASER_TIME = (LAYER_LENGTH + 2*H_SPACING) / LASER_SPEED   # Time on a pass")
    lines.append("LAYER_WIDTH = NB_LASER * H_SPACING  # Layer width (m)")
    lines.append("LAYER_TIME = REAL_LENGTH/LASER_SPEED - LASER_TIME # Time of cooling between each passes")
    lines.append("")

    # Material
    lines.append("# Material #")
    lines.append("")
    lines.append(f"MAT_NAME = '{material.get()}' # Material name")
    lines.append(f"LATENT_HEAT = {latent_heat}   # Latent heat of fusion (J/kg)")
    lines.append(f"SOLIDUS = {solidus}   # Solidus temperature (K)")
    lines.append(f"LIQUIDUS = {liquidus}   # Liquidus temperature (K)")

    lines.append(f"data_density = {data_density}   # Density (kg/m3)")
    lines.append(f"data_density_powder = {data_density_powder}   # Density (kg/m3)")
    lines.append(f"data_cond = {data_cond}   # Conductivity (W/m/K)")
    lines.append(f"data_spe_heat = {data_spe_heat}   # Specific Heat (J/kg/K)")
    lines.append(f"data_expansion = {data_expansion}   # Thermal Expansion (1/K)")
    lines.append(f"data_elastic = {data_elastic}   # Young Modulus (Pa), Poisson ratio")
    lines.append(f"data_plastic = {data_plastic}   # Yield strenght (Pa), Plastic Strain")
    lines.append(f"ambient_density = {ambient_density}   # Density at ambient temperature (kg/m3)")

    # Write file
    with open("constant.py", "w") as f:
        f.write("\n".join(lines))
    messagebox.showinfo("Success", "File saved:\nconstant.py")

# Write CONSTANT_UMATHT.inc based on input parameters and material tables
def write_umath_inc(porosity, solidus, liquidus):
    if data_cond is None or data_spe_heat is None:
        messagebox.showerror("Missing data", "Please import Material table first")
        return

    lines = []
    # Write scalar parameters
    lines.append(f"DOUBLE PRECISION, PARAMETER :: POROSITY = {porosity}")
    lines.append(f"DOUBLE PRECISION, PARAMETER :: SOLIDUS_TEMP = {solidus}")
    lines.append(f"DOUBLE PRECISION, PARAMETER :: LIQUIDUS_TEMP = {liquidus}")
    lines.append("")
    lines.append("! = Material properties = !")
    lines.append(f"DOUBLE PRECISION, DIMENSION(2, {len(data_cond)}) :: TABLE_K, TABLE_C")
    lines.append(f"DOUBLE PRECISION, DIMENSION(2, {len(data_cond_powder)}) :: TABLE_K_P ! Conductivity of powder")

    # Write conductivity table
    lines.append("")
    lines.append("! Thermal Conductivity")
    lines.append("DATA TABLE_K / &")
    for cond, temp in data_cond:
        lines.append(f" {float(cond)}, {float(temp)}, &") 
    lines[-1] = lines[-1].replace(', &', ' /')  # Remove trailing '&' from the last line
    lines.append("")
    lines.append("DATA TABLE_K_P / &")
    for cond, temp in data_cond_powder:
        lines.append(f" {float(cond)}, {float(temp)}, &") 
    lines[-1] = lines[-1].replace(', &', ' /')  # Remove trailing '&' from the last line

    # Write specific heat table
    lines.append("")
    lines.append("! Specific Heat")
    lines.append("DATA TABLE_C / &")
    for spe_heat, temp in data_spe_heat:
        lines.append(f" {float(spe_heat)}, {float(temp)}, &") 
    lines[-1] = lines[-1].replace(', &', ' /')  # Remove trailing '&' from the last line
    
    # Write file
    with open("CONSTANT_UMATHT.inc", "w") as f:
        f.write("\n".join(lines))
    messagebox.showinfo("Success", "File saved:\nCONSTANT_UMATHT.inc")


# Write CONSTANT_DFLUX.inc 
def write_dflux_inc(entry):
    lines = []
    for key in entry.keys():
        if key == 'NB_LASER' or key == 'NB_LAYER' :
            lines.append(f"INTEGER, PARAMETER :: {key} = {int(entry[key].get())}")
        else:
            lines.append(f"DOUBLE PRECISION, PARAMETER :: {key} = {float(entry[key].get()):e}".replace('e', 'D'))
    with open("CONSTANT_DFLUX.inc", "w") as f:
        f.write("\n".join(lines))
    messagebox.showinfo("Success", "File saved:\nCONSTANT_DFLUX.inc")

############################## GUI ###################################

# ------- Inputs ----
entries_dflux = {}
entries_python = {}

fields_layer = [
    ("Nb of layers", "NB_LAYERS", 4),
    ("Nb of laser passes per layer", "NB_LASER", 10),
    ("Hatch spacing (m)", "H_SPACING", 80e-6),
    ("Real layer length (m)", "REAL_LENGTH", 0.01),
    ("Sample length (m)", "LAYER_LENGTH", 800e-6),
    ("Layer thickness (m)", "LAYER_THICKNESS", 40e-6),
    ("Support thickness (m)", "SUPPORT_THICKNESS", 0.5e-3)]
fields_laser = [
    ("Laser Power (W)", "LASER_POWER", 190),
    ("Laser radius (m)", "LASER_RADIUS", 45e-6),
    ("Laser speed (m/s)", "LASER_SPEED", 800e-3),
    ("Cooling time btw layers (s)", "COOLING_TIME", 10)] # 1.5s is the time needed to print 1 layer of 1mm² at 0.8m/s
fields_environment = [
    ("Ambient temperature (K)", "AMBIENT_TEMP", 350),
    ("Convection coef (W/m2)", "FILMCOEFF", 15)]
fields_material = [
    ("Emissivity", "EMISSIVITY", 0.45),
    ("Absorptivity", "ABSORB", 0.35),
    ("Molar mass (kg/mol)", "M", 0.05514),
    ("Latent heat of vaporization (J/kg)", "DHV", 7.41e6)
    ]


tk_root = Tk()
tk_root.title("Simulation Parameters")

# === TOP ROW: JOB & ENVIRONMENT ===
frame_top = Frame(tk_root)
frame_top.pack(padx=10, pady=10)

# ---- JOB ----
frame_job = Frame(frame_top)
frame_job.grid(row=0, column=0, padx=20)
Label(frame_job, text="=== JOB PARAMETERS ===", font=("Helvetica", 10, "bold")).grid(row=0, column=0, columnspan=2, sticky="w", pady=5)

Label(frame_job, text="Thermal Job name").grid(row=1, column=0, sticky="w", pady=2)
job = Entry(frame_job, width=20)
job.insert(0, "Job-80h-45rad")
job.grid(row=1, column=1, pady=2)

Label(frame_job, text="Number of processors").grid(row=2, column=0, sticky="w", pady=2)
entry = Entry(frame_job, width=20)
entry.insert(0, str(6))
entry.grid(row=2, column=1, pady=2)
entries_python["NB_CPU"] = entry

# ---- ENVIRONMENT ----
frame_env = Frame(frame_top)
frame_env.grid(row=0, column=1, padx=20)
Label(frame_env, text="=== ENVIRONMENT ===", font=("Helvetica", 10, "bold")).grid(row=0, column=0, columnspan=2, sticky="w", pady=5)

for i, (label_text, var_name, default) in enumerate(fields_environment):
    Label(frame_env, text=label_text).grid(row=i+1, column=0, sticky="w", pady=2)
    entry = Entry(frame_env, width=20)
    entry.insert(0, str(default))
    entry.grid(row=i+1, column=1, pady=2)
    entries_python[var_name] = entry
    if i == 0 : #ambienttemp
        ambient_temp = float(entry.get())


# === BOTTOM ROW: LAYER & LASER ===
frame_bottom = Frame(tk_root)
frame_bottom.pack(padx=10, pady=10)

# ---- LAYER ----
frame_layer = Frame(frame_bottom)
frame_layer.grid(row=0, column=0, padx=20)
Label(frame_layer, text="=== LAYER ===", font=("Helvetica", 10, "bold")).grid(row=0, column=0, columnspan=2, sticky="w", pady=5)
Label(frame_layer, text="Layer width = h_spacing*nb_passes", font=('TkDefaultFont', 9, 'italic')).grid(row=1, column=0, columnspan=2, sticky="w", pady=5)

for i, (label_text, var_name, default) in enumerate(fields_layer):
    Label(frame_layer, text=label_text).grid(row=i+2, column=0, sticky="w", pady=2)
    entry = Entry(frame_layer, width=20)
    entry.insert(0, str(default))
    entry.grid(row=i+2, column=1, pady=2)
    entries_dflux[var_name] = entry

# ---- LASER ----
frame_laser = Frame(frame_bottom)
frame_laser.grid(row=0, column=1, padx=20)
Label(frame_laser, text="=== LASER ===", font=("Helvetica", 10, "bold")).grid(row=0, column=0, columnspan=2, sticky="w", pady=5)

for i, (label_text, var_name, default) in enumerate(fields_laser):
    Label(frame_laser, text=label_text).grid(row=i+1, column=0, sticky="w", pady=2)
    entry = Entry(frame_laser, width=20)
    entry.insert(0, str(default))
    entry.grid(row=i+1, column=1, pady=2)
    entries_dflux[var_name] = entry


# ---- MATERIAL ----
Label(tk_root, text="=== MATERIAL ===", font=("Helvetica", 10, "bold")).pack(pady=5)
frame_mat = Frame(tk_root)
frame_mat.pack(padx=10, pady=10)

# Material Name #
Label(frame_mat, text="Material name").grid(row=0, column=0, sticky="w", pady=2)
material = Entry(frame_mat, width=20)
material.insert(0, "SS-316L")
material.grid(row=0, column=1, pady=2)

# Material properties #
Label(frame_mat, text="Porosity").grid(row=1, column=0, sticky="w", pady=2)
porosity = Entry(frame_mat, width=20)
porosity.insert(0, "0.4")
porosity.grid(row=1, column=1, pady=2)

Label(frame_mat, text="Solidus temperature (K)").grid(row=2, column=0, sticky="w", pady=2)
solidus = Entry(frame_mat, width=20)
solidus.insert(0, "1643.0")
solidus.grid(row=2, column=1, pady=2)

Label(frame_mat, text="Liquidus temperature (K)").grid(row=3, column=0, sticky="w", pady=2)
liquidus = Entry(frame_mat, width=20)
liquidus.insert(0, "1700.0")
liquidus.grid(row=3, column=1, pady=2)

Label(frame_mat, text="Latent heat of fusion (J/kg)").grid(row=4, column=0, sticky="w", pady=2)
latent_heat = Entry(frame_mat, width=20)
latent_heat.insert(0, "290000")
latent_heat.grid(row=4, column=1, pady=2)

for i, (label_text, var_name, default) in enumerate(fields_material):
    Label(frame_mat, text=label_text).grid(row=5+i, column=0, sticky="w", pady=2)
    entry = Entry(frame_mat, width=20)
    entry.insert(0, str(default))
    entry.grid(row=5+i, column=1, pady=2)
    entries_dflux[var_name] = entry

# Material Table #
def import_therm_table():
    global data_therm 
    data_therm = import_csv()
    if data_therm:
        read_therm_table(data_therm, float(porosity.get()), float(solidus.get()), float(liquidus.get()), ambient_temp)
        
def import_mecha_table():
    global data_mecha 
    data_mecha = import_csv()
    if data_mecha:
        read_mecha_table(data_mecha)
        
def import_pla_table():
    global data_pla 
    data_pla = import_csv()
    if data_pla:
        read_pla_table(data_pla)

frame_import = Frame(frame_mat)
frame_import.grid(row=len(fields_material)+5, column=0, columnspan=2, pady=5)

# Ligne des boutons
frame_buttons = Frame(tk_root)
frame_buttons.pack(pady=10)
Button(frame_buttons, text="Import Thermal table", command=import_therm_table).pack(side='left', padx=10)
Button(frame_buttons, text="Import Mechanic table", command=import_mecha_table).pack(side='left', padx=10)
Button(frame_buttons, text="Import Plastic table", command=import_pla_table).pack(side='left', padx=10)

# Label en dessous
Label(frame_import,text="Complete thermal, mechanic & plastic.csv with material properties values",fg='red', font=('TkDefaultFont', 9, 'italic')).pack(pady=5)

# ---- Save buttons ----
frame_btn = Frame(tk_root)
frame_btn.pack(pady=10)
Button(frame_btn, text="💾 Save CONSTANT_DFLUX.inc", command=lambda: write_dflux_inc(entries_dflux), bg="black", fg="white").pack(side='left', padx=10)
Button(frame_btn, text="💾 Save CONSTANT_UMATHT.inc", command=lambda: write_umath_inc(float(porosity.get()), float(solidus.get()), float(liquidus.get())), bg="black", fg="white").pack(side='left', padx=10)
Button(frame_btn, text="💾 Save constant.py", command=lambda: write_constant_py(job, material, float(latent_heat.get()), float(solidus.get()), float(liquidus.get()), entries_dflux|entries_python), bg="black", fg="white").pack(side='left', padx=10)

# Start the GUI
tk_root.mainloop()
