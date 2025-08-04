
#files 
import Surface_rays as geo
import Surface_data as FuD
import NBI_Ports_data_input as Cout
import J_0_test.mconf.mconf as mconf
import Weight_Fuction.WF_FIDA_CTS as WF
import subprocess
import sys
import platform
from scipy.integrate import trapezoid

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

libraries = [ "tkinter", "matplotlib","customtkinter", "numpy","scipy", "concurrent.futures","tqdm", "json", "jsonpickle"]
for library in libraries:
    try:
        globals()[library] = __import__(library)  
    except ImportError:
        print(f"Loading....: {library}")
        install(library)  


#library 
import tkinter as tk
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import customtkinter as ctk
import numpy as np
import os 
from scipy.integrate import solve_ivp
from matplotlib import cm
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import json
from tkinter import filedialog
import jsonpickle
from concurrent.futures import ThreadPoolExecutor




class Data:
        #====================================================================DATA TYPE ====================================================================================================
        #data_B: [Name NBI; Name Port; Points on NBI; Mag field in this points; Angle between linesight and vec NBI; vec Mag field in points on NBI; angle between vec linesi and magfield]
        #data_B[0]: Name NBI
        #data_B[1]: Name Port
        #data_B[2]: Points on NBI
        #data_B[3]: Mag field in this points
        #data_B[4]: Angle between linesight and vec NBI
        #data_B[5]: vec Mag field in points on NBI
        #data_B[6]: angle between vec linesi and magfield
        #data_B[7]: S
        #data_B[8]: B_max 
        #data_B[9]: J_0
        #data_B[10]: WF
        #==================================================================================================================================================================================



    def __init__(self):
        self.R_x, self.R_y, self.R_z = FuD.all_point(FuD.read_data()[0])
        self.P_1, self.P_2, self.P_name = Cout.Ports()
        self.NBI_X, self.NBI_Y, self.NBI_Z, self.NBI_uvec_X, self.NBI_uvec_Y, self.NBI_uvec_Z = Cout.NBI()
        self.Bget = calculus()
        self.resolution_WF = 300



        # ===================================================================================================Create 3D surface=========================================================
        self.surface = geo.create_surface(self.R_x, self.R_y, self.R_z)

        # Get intersections for ports and NBI
        self.new_P_1, *_ = geo.get_intersection_points(self.P_1, self.P_2, self.surface)
        self.new_NBI_start, self.new_NBI_end, *_ = geo.get_intersection_points_NBI(
            self.NBI_X, self.NBI_Y, self.NBI_Z, self.NBI_uvec_X, self.NBI_uvec_Y, self.NBI_uvec_Z, self.surface
        )
        self.valid_indx = []
        for i in range(0, 11):
            valid_indices = geo.pre_NBI_and_PORTS(i, self.new_P_1, self.new_NBI_start, self.new_NBI_end, self.surface)
            # Store the valid port names
            self.valid_indx.append(valid_indices)
        #============================================================================================================================================================================

    def port_for_nbi(self, NBI_index, angle, scale):
        P_1_for_NBI = self.new_P_1[:, self.valid_indx[NBI_index]]
        P_1_start_for_NBI = self.P_1[:, self.valid_indx[NBI_index]]
        Pname_for_NBI = [self.P_name[i] for i in self.valid_indx[NBI_index]]
        valid_indices, extreme_points_1, extreme_points_2, *_ = geo.NBI_and_PORTS(
            P_1_start_for_NBI, NBI_index, P_1_for_NBI, self.new_NBI_start, self.new_NBI_end, self.surface, float(angle))
        valid_port_names = [Pname_for_NBI[i] for i in valid_indices]
        return valid_indices, extreme_points_1, extreme_points_2, valid_port_names
    

    def data_already_input(self, scale, Name_Ports, Name_NBI, angle, config, B_0, diagnostic_types):
      
        data_B = [[] for _ in range(11)]
        for i in range(len(Name_Ports)):
           data_i = self.data_nbi_ports(Name_NBI[i], Name_Ports[i], angle, scale, config,B_0)
           
           for i in range(len(data_i)):
               data_B[i].append(data_i[i])
           print(data_B[0])

        WF_array = self.create_result_array_for_port_old(data_B, diagnostic_types)
        data_B[10].append(WF_array)

        return data_B
    


    def data_nbi_ports(self, nbi, port, angle, scale, config, B_0):
        index = self.P_name.index(port)
        P_1_start = [self.P_1[0][index], self.P_1[1][index], self.P_1[2][index]]
        P_2_end = [self.new_P_1[0][index], self.new_P_1[1][index], self.new_P_1[2][index]]
        
        index_NBI = int(nbi.split('_')[-1])-1
        if nbi.startswith("Gyrotron"):
            index_NBI += 8

        valid_indices, extreme_points_1, extreme_points_2, *_ = geo.NBI_and_PORTS(
            P_1_start, index_NBI, P_2_end, self.new_NBI_start, self.new_NBI_end, self.surface, float(angle))

        points, B_array, B_vec_array, S_array, B_max_array= self.Bget.gets(np.array(extreme_points_1[0], dtype=np.float64), np.array(extreme_points_2[0], dtype=np.float64), scale, config, B_0)

        angles, angles_vec_B, J_0_array=[],[], []
        for j in range(len(points)):
                angle = geo.check_segment_angle(P_1_start, P_2_end, points[j]*100)
                angles.append(angle)
                vector_AB = np.array(P_2_end) - np.array(P_1_start)  
                angle_B = geo.check_angle_2_vec(vector_AB/100, B_vec_array[j])
                angles_vec_B.append(angle_B)
        J_0_array = self.Bget.J_0_calculate(points, config, B_0)

        return [nbi, port, points, B_array, angles, B_vec_array, angles_vec_B,S_array, B_max_array, J_0_array] 


    def create_result_array_for_port_old(self, data_B_c, diagnostic_types):
        WF_array = []
        import pickle
        pickle.dumps(data_B_c) 
        with ProcessPoolExecutor(max_workers=5) as executor:
            results = list(tqdm(executor.map( self.Bget.process_data_for_point, range(len(data_B_c[1])), [data_B_c] * len(data_B_c[0]), diagnostic_types)))

        for result_for_i in results:
            WF_array.append(result_for_i)
        print("ready!")
    
        return results
""""
    def process_data_for_point(self,i, data_B_c, diagnostic_type):

        index_nbi = int(data_B_c[0][i].split('_')[-1]) - 1
        result_for_j = []

        for j in range(len(data_B_c[3][i])):
            x_ev = np.linspace(10, 100, self.resolution_WF)
            y_ev = np.linspace(-100, 100, self.resolution_WF) / 2.5

            if diagnostic_type == 'FIDA':
               result = WF.weight_Function(data_B_c[4][i][j], data_B_c[3][i][j], x_ev, y_ev)
            if diagnostic_type == 'CTS':
                result = WF.CTS_wf(data_B_c[6][i][j], data_B_c[3][i][j], x_ev, y_ev)
            result_for_j.append(result) 
        return result_for_j
"""



class calculus():
    def __init__(self):
         lll=0
         self.resolution_WF = 300



    def process_data_for_point(self, i, data_B_c, diagnostic_type):
        index_nbi = int(data_B_c[0][i].split('_')[-1]) - 1
        result_for_j = []

        for j in range(len(data_B_c[3][i])):
            x_ev = np.linspace(10, 100, self.resolution_WF)
            y_ev = np.linspace(-100, 100, self.resolution_WF) / 2.5

            if diagnostic_type == 'FIDA':
                result = WF.weight_Function(data_B_c[4][i][j], data_B_c[3][i][j], x_ev, y_ev)
            elif diagnostic_type == 'CTS':
                result = WF.CTS_wf(data_B_c[6][i][j], data_B_c[3][i][j], x_ev, y_ev)
            result_for_j.append(result)
        return result_for_j
    


    def gets(self, point1, point2, scale, config, B_0):
      point1, point2 = point1 / 100, point2 / 100
      points = np.linspace(point1, point2, scale)

      previous_directory = os.getcwd()
      print(previous_directory)
      os.chdir('J_0_test')


      mconf_config = {'B0': B_0,
                'B0_angle': 0.0,
                'accuracy': 1e-10, 
                'truncation': 1e-10} 
      eq = mconf.Mconf_equilibrium(config ,mconf_config=mconf_config)

      
      S_array = []
      for i in range(len(points)):
         S, vecB = eq.get_B(points[i])
         S_array.append(S)
      start_indices, end_indices = self.find_transitions(S_array)
      start_point = points[start_indices[0]].reshape(3,)
      print(start_indices)
      print(end_indices)
      end_point = points[end_indices[0]].reshape(3,)

      points = np.linspace(start_point, end_point, scale)

      B_array, B_vec_array, S_array, B_max_array= [], [], [], []
      for i in range(len(points)):
         S, vecB = eq.get_B(points[i])
         B_max = eq.get_Bmax(S)
         valueB = np.sqrt(vecB[0]**2 + vecB[1]**2 + vecB[2]**2)
         S_array.append(S)
         B_array.append(valueB)
         B_vec_array.append(vecB)
         S_array.append(S)
         B_max_array.append(B_max)

      os.chdir(previous_directory)
      return points, B_array, B_vec_array, S_array, B_max_array
    

    def find_transitions(self, arr):
        start_indices = []
        end_indices = []
    
        for i in range(len(arr) - 1):
            if arr[i] >= 1 and arr[i+1] < 1:   
                start_indices.append(i+1)
            elif arr[i] < 1 and arr[i+1] >= 1: 
                end_indices.append(i+1)
        if not start_indices:
                start_indices.append(0)
        if not end_indices:
                end_indices.append(len(arr) - 1)

    
        return start_indices, end_indices







    #------------------------------------------------------------------------------J_0 calc-------------------------------------------------------------------------------------------------
    def J_0_calculate(self, points,config, B_0):
     points = np.array(points, dtype=np.float64)
     previous_directory = os.getcwd()
     os.chdir('J_0_test')
    
     with ProcessPoolExecutor(max_workers=5) as executor:
        results = list(tqdm(executor.map(self.calculate_J_0_for_point, points, [config] * len(points), [B_0]*len(points)), total=len(points)))
     print(config)
     os.chdir(previous_directory)
     return results

    def calculate_J_0_for_point(self, point, config, B_0):
      #data type
      point = np.array(point, dtype=np.float64)

      #config
      mconf_config = {'B0': B_0,
                'B0_angle': 0.0,
                'accuracy': 1e-10, #accuracy of magnetic to cartesian coordinat transformation
                'truncation': 1e-10} #trancation of mn harmonics
      

      eq = mconf.Mconf_equilibrium(config,mconf_config=mconf_config)
      

      #constant
      s0, vecB = eq.get_B(point)
      B_value = np.linalg.norm(vecB)
      B_max_point = eq.get_Bmax(s0)
      L = 300 
      N = 2000
      E_values = np.linspace(10, 100, self.resolution_WF) * (1.6 * 10**(-19)) * 10**3  
      mu_values = np.linspace(-100, 100, self.resolution_WF)/2.5 * (1.6 * 10**(-19)) * 10**3 

      #Solve eq
      def solve_differential(point, L, N, rhs_B):
        sol = solve_ivp(rhs_B, [0, L], point,method='RK45', max_step=L / N,atol=1e-6, dense_output=True)
        s_sol, B_sol = eq.get_s_B_T(sol.y[0],sol.y[1], sol.y[2])
        magB = np.linalg.norm(B_sol, axis=1)
        path = np.zeros(sol.y.shape[1])
        path = np.cumsum(np.sqrt(np.sum(np.diff(sol.y,axis=-1)**2, axis=0)))
        path = np.insert(path, 0, 0) 
        return magB, path

      
      rhs_B_forward = lambda l, y: eq.get_B(y)[1]/ np.linalg.norm(eq.get_B(y)[1])
      rhs_B_backward = lambda l, y: -eq.get_B(y)[1]/ np.linalg.norm(eq.get_B(y)[1])

      forward_magB, forward_path = solve_differential(point, L, N, rhs_B_forward)
      backward_magB, backward_path = solve_differential(point, L, N, rhs_B_backward)

    
      #Arrays
      E_values, mu_values = np.meshgrid(E_values, mu_values)
      J_0_map = np.zeros(E_values.shape)

      for i in range(E_values.shape[0]):
          for j in range(mu_values.shape[1]):
              
              E = E_values[i, j]
              mu = mu_values[i, j]
              B_max_particle = np.abs(E / mu)
              
              
              if B_max_particle<=B_max_point and B_max_particle>B_value and s0<=1:
               
               #-------------------Integral----------------------------
               def compute_integrals(magB, path, mu, B_max_particle):
                 mask = magB <= B_max_particle
                 idx_limit = np.argmax(~mask) if np.any(~mask) else len(magB)
                 magB_limited, path_limited = magB[:idx_limit], path[:idx_limit]
                 integrand = np.sqrt(2 * np.abs(mu) * (B_max_particle - magB_limited))
                 return integrand, path_limited
    

               forward_integrand, forward_path_limited = compute_integrals(forward_magB, forward_path, mu, B_max_particle)
               backward_integrand, backward_path_limited = compute_integrals(backward_magB, backward_path, mu, B_max_particle)
              

               complete_path = np.concatenate([
                backward_path_limited[::-1],  
                forward_path_limited,         
                forward_path_limited[::-1],   
                backward_path_limited         
                ])

               complete_integrand = np.concatenate([
                 backward_integrand[::-1],    
                 forward_integrand,           
                 forward_integrand[::-1],     
                 backward_integrand          
                 ])
               
               J_0_map[i, j] = self.trapezoidal_integral(complete_integrand, complete_path)*1e7
              else:
               J_0_map[i, j] = np.nan

      return J_0_map
    
    def trapezoidal_integral(self, f, s):
      ds = np.abs(np.diff(s))
      avg_f = (f[:-1] + f[1:]) / 2
      segment_integrals = ds * avg_f
      integral = np.sum(segment_integrals)
      return integral
    


    def compute_matrix(self, args, all_results, delta_J, delta_s):
        i, j, port_i, port_j = args
        matrix = self.sum(port_i, port_j, i, j, all_results, delta_J, delta_s)
        #matrix = np.transpose(matrix)
        filtered = matrix[matrix != -np.inf]
        min_value = np.min(filtered) if filtered.size > 0 else 0
        return (i, j, matrix, min_value)




    def sum(self, array_1, array_2, first_nbi_index, second_nbi_index, all_results, delta_J, delta_s):
     MATRIX = np.zeros((len(array_1), len(array_2)))

     
     x_ev = np.linspace(10, 100, self.resolution_WF)
     y_ev = np.linspace(-100, 100, self.resolution_WF) / 2.5
     x_ev, y_ev = np.meshgrid(x_ev, y_ev)

     ratio = np.abs(x_ev / y_ev)

     B_value_1 = all_results[3][first_nbi_index]
     B_value_2 = all_results[3][second_nbi_index]
     B_max_array_1 = all_results[8][first_nbi_index]
     B_max_array_2 = all_results[8][second_nbi_index]
     s_1_array = all_results[7][first_nbi_index]
     s_2_array = all_results[7][second_nbi_index]
     J_0_array_1 = np.array(all_results[9][first_nbi_index])  
     J_0_array_2 = np.array(all_results[9][second_nbi_index])  
    
     for i in range(len(array_1)):
        for j in range(len(array_2)):

            B_max_i = B_max_array_1[i]
            B_max_j = B_max_array_2[j]
            s_1_i = s_1_array[i]
            s_2_j = s_2_array[j]

            mask_i = (ratio > B_max_i)
            mask_j = (ratio > B_max_j)

            mask_i_val = (ratio> B_value_1[i])
            mask_j_val = (ratio> B_value_2[j])
            
            #free particles 
            both_above_B_mask = mask_i & mask_j

            #not free
            both_below_B_mask = np.logical_not(mask_i) & np.logical_not(mask_j) & mask_i_val & mask_j_val

            #or or 
            cross_check_mask = (mask_i & np.logical_not(mask_j)) | (np.logical_not(mask_i) & mask_j)


            
            
            
            
            #==================================NOT free=================================================
            def relative_difference(a, b):
                numerator = np.abs(a - b)
            
                denominator = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-13)
                return numerator / denominator
            nan_mask = np.logical_or(np.isnan(J_0_array_1[i]), np.isnan(J_0_array_2[j])) 
            equal_mask_J_0 = np.where(nan_mask, True, relative_difference(J_0_array_1[i], J_0_array_2[j])<delta_J)
            
            
            #==================================free=================================================
            equal_mask_s = (np.abs(np.abs(s_1_i / s_2_j)-1) <=  delta_s)
            

            #==================================forbidden======================================
            mask_no_accept =  np.logical_not(mask_i_val) | np.logical_not(mask_j_val) 

            #free diff
            mask_condition_1 = both_above_B_mask  & np.logical_not(equal_mask_s)
            #free same
            mask_condition_2 = both_above_B_mask & equal_mask_s
            #or or 
            mask_condition_4 = cross_check_mask
            #not free same 
            mask_condition_5 = both_below_B_mask & equal_mask_J_0
            #not free diff 
            mask_condition_6 = both_below_B_mask & np.logical_not(equal_mask_J_0)
            #f
            mask_condition_7 = mask_no_accept
            
            #results
            product = array_1[i] * array_2[j]
            #mask
            product[mask_condition_1] = 0
            product[mask_condition_4] = 0
            product[mask_condition_6] = 0
            product[mask_condition_7] = 0




            sum_product = np.trapz(np.trapz(product, np.linspace(10, 100, self.resolution_WF), axis=1), np.linspace(-100, 100, self.resolution_WF) / 2.5)
            element = np.log10(np.where(sum_product > 0, sum_product, 1e-8))
            MATRIX[i, j] = np.where(element>-8, element, -np.inf)

     return MATRIX

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
