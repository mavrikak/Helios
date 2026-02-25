import numpy as np
import os
import sys

# --- Physical Constants for Unit Conversion ---
HC_EV_NM = 1239.84193

def calculate_drude_dielectric(wavelengths_nm, eps_inf, E_p_eV, E_gamma_eV):
    # Convert wavelength (nm) to Energy (eV)
    with np.errstate(divide='ignore'):
        E = HC_EV_NM / wavelengths_nm
        E[np.isinf(E)] = 0 

    E_p_sq = E_p_eV**2
    E_gamma_sq = E_gamma_eV**2
    Denominator = E**2 + E_gamma_sq

    e1 = eps_inf - (E_p_sq / Denominator)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        e2 = (E_p_sq * E_gamma_eV) / (E * Denominator)
        e2 = np.nan_to_num(e2)

    return e1, e2

def save_material_file(filename, material_title, wavelengths, e1, e2):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.abspath(os.path.join(script_dir, os.pardir, 'materials'))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    full_path = os.path.join(output_dir, filename)
    print("\nSaving file to: {}".format(full_path))

    try:
        with open(full_path, 'w') as f:
            header = "real and imaginary epsilon of {} as a function of wavelength (nm), Drude Model Calculation\n".format(material_title)
            f.write(header)
            
            for w, re, im in zip(wavelengths, e1, e2):
                line = "{:.2f}\t{:.7f}\t{:.7f}\n".format(w, re, im)
                f.write(line)
        print("File saved successfully.")
        
    except IOError as e:
        print("Error writing file: {}".format(e))

def main():
    print("--- Drude Model Material Generator ---\n")
    print("This script generates optical constants using the Drude model.")
    
    if sys.version_info[0] < 3:
        input_func = raw_input
    else:
        input_func = input

    try:
        material_name = input_func("Enter material name for header (e.g., 'Gold'): ").strip()
        filename = input_func("Enter desired filename (e.g., 'Au_Drude.txt'): ").strip()
        if not filename.endswith('.txt'):
            filename += ".txt"

        print("\n-- Drude Parameters --")
        eps_inf = float(input_func("Enter Epsilon Infinity (epsilon_inf): "))
        E_p = float(input_func("Enter Plasma Energy (E_p) in eV: "))
        E_gamma = float(input_func("Enter Damping Energy (Gamma) in eV: "))

        print("\n-- Wavelength Range (nm) --")
        w_start = float(input_func("Enter start wavelength (nm): "))
        w_end = float(input_func("Enter end wavelength (nm): "))
        num_points = int(input_func("Enter number of data points: "))

        if w_start <= 0 or w_end <= 0:
             print("Error: Wavelengths must be positive.")
             return

    except ValueError:
        print("\nError: Invalid numerical input.")
        return

    print("\nCalculating...")
    wavelengths = np.linspace(w_start, w_end, num_points)
    e1, e2 = calculate_drude_dielectric(wavelengths, eps_inf, E_p, E_gamma)

    save_material_file(filename, material_name, wavelengths, e1, e2)

if __name__ == "__main__":
    main()