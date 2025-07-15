import numpy as np
import pandas as pd
import os


def NBI():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_dir, 'Input_data', 'NBI_Coordinates.xlsx')
    df = pd.read_excel(file_path)

    #NBI_ports_start
    NBI_X = df.iloc[4:, 4].tolist()  # 5 colum, from 6 row
    NBI_Y = df.iloc[4:, 5].tolist()  # 6 colum, from 6 row
    NBI_Z = df.iloc[4:, 6].tolist()  # 7 colum, from 6 row
    
    #Vector
    NBI_uvec_X = df.iloc[4:, 7].tolist()  # 8 colum, from 6 row
    NBI_uvec_Y = df.iloc[4:, 8].tolist()  # 9 colum, from 6 row
    NBI_uvec_Z = df.iloc[4:, 9].tolist()  # 10 colum, from 6 row
    
    
    #float
    NBI_X = [float(value) for value in NBI_X]
    NBI_Y = [float(value) for value in NBI_Y]
    NBI_Z = [float(value) for value in NBI_Z]
    
    #float
    NBI_uvec_X = [float(value) for value in NBI_uvec_X]
    NBI_uvec_Y = [float(value) for value in NBI_uvec_Y]
    NBI_uvec_Z = [float(value) for value in NBI_uvec_Z]
    

    
    return NBI_X, NBI_Y, NBI_Z, NBI_uvec_X, NBI_uvec_Y, NBI_uvec_Z


def Ports():
        current_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(current_dir, 'PortsDATA_ver2.xlsx')
        df = pd.read_excel(file_path)
    # P_1

        P_1_X = df.iloc[2:, 4].tolist() #4 row 5 colum 
        P_1_Y = df.iloc[2:, 5].tolist() #4 row 6 colum 
        P_1_Z = df.iloc[2:, 6].tolist() #4 row 7 colum
        
        P_1_X = [float(value) for value in P_1_X]
        P_1_Y = [float(value) for value in P_1_Y]
        P_1_Z = [float(value) for value in P_1_Z]
        

        
        P_1 = [P_1_X, P_1_Y, P_1_Z]
        P_1 = np.array(P_1)

    # P_2

        P_2_X = df.iloc[2:, 7].tolist() #4 row 8 colum 
        P_2_Y = df.iloc[2:, 8].tolist() #4 row 9 colum 
        P_2_Z = df.iloc[2:, 9].tolist() #4 row 10 colum
        
        P_2_X = [float(value) for value in P_2_X]
        P_2_Y = [float(value) for value in P_2_Y]
        P_2_Z = [float(value) for value in P_2_Z]
        
        
        P_2 = [P_2_X, P_2_Y, P_2_Z]
        P_2 = np.array(P_2)




    #ports_name
        P_name = df.iloc[2:, 1].tolist()  
        P_module = df.iloc[2:, 2].tolist()  
        P_submodule = df.iloc[2:, 3].tolist()  

        P_name = np.array(P_name)
        P_module = [int(value) for value in P_module]

        # Remove the first character from each port name
        modified_P_name = [name[2:] for name in P_name]  # Keep only the characters after the first one

        # Convert P_submodule to the desired format
        modified_P_submodule = [0 if submodule == 0 else 1 for submodule in P_submodule]

        # Create the new array with the desired format
        combined_array = [f"{module}_{submodule}_{name}" for module, submodule, name in zip(P_module, modified_P_submodule, modified_P_name)]


        return P_1, P_2, combined_array


