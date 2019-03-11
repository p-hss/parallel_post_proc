import math
def read_data(i,j,k):

    sim=["R01", "R02", "R03", "R04",
         "R05", "R06", "R07", "R08", ]

    files_1=["../data/sim_",
             "/home/gamling/p.hess/Masterarbeit/Thesis"]
    files_2=["/tables/tab_results_", 
            "/tables/tab_results_angle_line_",
            "/tables/tab_mhd_sim_charact_"]
    files_3=[".txt", "_1.txt"]
    if k == 1:
        file=files_1[k]+files_2[j]+sim[i]+files_3[k]
    else:
        file=files_1[k]+sim[i]+files_2[j]+sim[i]+files_3[k]

    f = open(file, "r")
    data=f.read()

    return data


file_out="../data/line_xhel_scaling.dat"
start=7
stop=0

for sims in range(start,stop,-1):

    data=read_data(sims,2,1)
    sigma=math.cos(float(data[61:66]))

    data=read_data(sims,0,0)
    zeta=data[40:45]
    zeta_error=data[46:54]
    if sims==start:
        max_zeta=float(zeta)

    data=read_data(sims,1,0)
    angle_max_strain=data[9:14]
    angle_max_strain_error=data[15:23]
    if sims==start:
        max_strain_angle=float(angle_max_strain)

    data=read_data(sims,1,0)
    angle_b=data[89:94]
    angle_b_error=data[95:103]
    if sims==start:
        max_b_angle=float(angle_b)

    output=str(sigma)+" "+str(float(zeta)/max_zeta)+" "+zeta_error+" "+ \
                          str(float(angle_max_strain)/max_strain_angle)+" "+angle_max_strain_error+" "+ \
                          str(float(angle_b)/max_b_angle)+" "+angle_b_error+" \n"
    
    if sims == start:
        f = open(file_out, "w")
        f.write(output)

    else: 
        f = open(file_out, "a")
        f.write(output)

    print(output)