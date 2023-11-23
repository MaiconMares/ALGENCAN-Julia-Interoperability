import glob, re, os

ALGENCAN_PATH = os.environ['HOME'] + '/Desktop/ALGENCAN-Julia-Interoperability'
problems_filename = ALGENCAN_PATH + "/problems.txt"
problems_file = open(problems_filename,'r')
problems_content = problems_file.read()
problems_list = re.split("\.SIF\s*",problems_content)[:-1]
problems_file.close()

for problem in problems_list:
    os.chdir(ALGENCAN_PATH + '/algencan-4.0.0')
    make_cmd = "make cutest PROBNAME=" + problem
    os.system(make_cmd)

    os.chdir(ALGENCAN_PATH + "/algencan-4.0.0/bin/cutest")
    execute_binary_cmd = "./algencan"
    os.system(execute_binary_cmd)
    
    result_filename = 'tabline.txt'
    curr_file = open(result_filename,'r')
    problem_result = curr_file.read()
    curr_file.close()
    
    filename = 'results.txt'
    file = open(filename,'a')

    content = "{}\t\t\t\t{}".format(problem, problem_result)
    file.write(content)
    file.close()