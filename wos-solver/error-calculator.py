import numpy, sys, math, csv


def read_csv_files(path):
    try:
        with open(path, 'r', newline='') as file:
            points = csv.reader(file)
            pointData = list(points)
            pointData = convert_str_to_float(pointData)
        return pointData
    except FileNotFoundError:
        print(f"Could not find file {path}")
        return None, None
    except Exception as e:
        print(f"Unknown error {e}")
        return None, None

def convert_str_to_float(pointData):
    if not pointData:
        return []
    
    converted_data = []
    for row in pointData:
        converted_row = []
        for item in row:
            try:
                val_double = float(item)
                converted_row.append(val_double)
            except ValueError as e:
                print(f"Could not convert data {item}. Continuing...")
                return None
        converted_data.append(converted_row)
    return converted_data    

def calculate_l1_error(analytical, wos):
    l1_error = 0
    for row in range(len(analytical)):
        _, _, analytical_T = analytical[row]
        wos_T = wos[row][2] # ASSUME they have same x,y, rows
        l1_error += abs(analytical_T - wos_T)
    
    l1_error = l1_error / len(analytical) # error = 1/n * sum( abs(error))
    return l1_error

def calculate_l2_error(analytical, wos):
    l2_error = 0
    for row in range(len(analytical)):
        _, _, analytical_T = analytical[row]
        wos_T = wos[row][2] # ASSUME they have same x,y, rows
        l2_error += abs(analytical_T - wos_T)**2
    
    l2_error = math.sqrt(l2_error / len(analytical)) # error = 1/n * sum( abs(error))
    return l2_error

def main():
    # read input arguments from cmd line
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <analytical_solution_path.csv> <wos_solution_path.csv>")
        sys.exit(1)
    
    analytical_path = sys.argv[1]
    wos_path = sys.argv[2]

    analytical_results = read_csv_files(analytical_path)
    wos_results = read_csv_files(wos_path)

    # both are valid, perform L1 and L2 comparison
    if analytical_results and wos_results and len(analytical_results) == len(wos_results):
        # L1 - Mean Absolute Error
        l1_error = calculate_l1_error(analytical_results, wos_results)
        
        # L2 - Root Means Squared Error
        l2_error = calculate_l2_error(analytical_results, wos_results)
        
        print(f"L1 error = {l1_error}")
        print(f"L2 error = {l2_error}")
    
if __name__ == "__main__":
    main()