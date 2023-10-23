import glob, re, os

ALGENCAN_PATH = os.environ['HOME'] + '/Desktop/ALGENCAN-Julia-Interoperability'
files_path = "../algencan-4.0.0/sources/interfaces/cutest/*.SIF"
found_problem_paths = glob.glob(files_path)

for file in found_problem_paths:
    splitted_dirs = re.split("\/", file)
    problem_name = re.split("\.SIF", splitted_dirs[-1])[0]

    os.chdir(ALGENCAN_PATH + '/algencan-4.0.0')
    make_cmd = "make cutest PROBNAME=" + problem_name
    os.system(make_cmd)

    os.chdir(ALGENCAN_PATH + "/algencan-4.0.0/bin/cutest")
    execute_binary_cmd = "./algencan"
    os.system(execute_binary_cmd)
    
    new_file_name = problem_name + '_tabline.txt'
    os.rename('tabline.txt',new_file_name)

    file = open(new_file_name,'r')
    lines = file.readlines()

    file.close()
    cpu_time = lines[0].split()[2]
    computed_obj_func = -1
    
    idx = 3
    for value in lines[0].split()[3:-1]:
        if 'D' in value:
            computed_obj_func = value
            break
        
        idx += 1
    
    filename = 'results.txt'
    file = open(filename,'a')

    content = "{} {} {}\n".format(problem_name, cpu_time, computed_obj_func)
    file.write(content)
    file.close()