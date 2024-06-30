import numpy as np
import math
from fractions import Fraction as fr

def gomory_cut_algo():
    input_file = open(r"C:\Users\naman\OneDrive\Desktop\College\Sem 4\MTL103\test_case_ass_2\test_case_10.txt", "r")
    data = input_file.readlines()
    precision = 0

    boo = 1 if (data[1].strip().lower() == "maximize") else 0
    data[4] = data[4][:-1].split()
    no_of_vars = len(data[4])
    for j in range(no_of_vars):
        if data[4][j][-1] == ",":
            data[4][j] = data[4][j][:-1]
        data[4][j] = fr(data[4][j])
    a = np.array(data[4], dtype = fr).reshape((1, -1))

    i = 5
    while (data[i] != '\n'):
        data[i] = data[i][:-1].split()
        for j in range(no_of_vars):
            if data[i][j][-1] == ",":
                data[i][j] = data[i][j][:-1]
            data[i][j] = fr(data[i][j])
        temp = np.array(data[i], dtype = fr).reshape((1, -1))
        a = np.concatenate((a, temp), axis = 0, dtype = fr)
        i += 1
    m, n = a.shape
    i += 2

    b = np.full((m, 1), fr("0"), dtype = fr)
    for j in range(m):
        b[j, 0] = fr(data[i])
        i += 1
    i += 2

    aug_lst = []
    cnt = 0
    for j in range(m):
        if data[i].strip() == ">=":
            a[j] *= fr("-1")
            b[j] *= fr("-1")
            aug_lst.append(cnt)
            cnt += 1
        elif data[i].strip() == "<=":
            aug_lst.append(cnt)
            cnt += 1
        else:
            aug_lst.append(-1)
        i += 1
    aug_mat = np.full((m, cnt), fr("0"), dtype = fr)
    for j in range(m):
        if aug_lst[j] != -1:
            aug_mat[j, aug_lst[j]] = fr("1")
    a = np.concatenate((a, aug_mat), axis = 1, dtype = fr)
    i += 2

    for j in range(m):
        if b[j, 0] < 0:
            a[j, :] *= fr("-1")
            b[j, 0] *= fr("-1")

    data[i] = data[i].strip().split()
    for j in range(no_of_vars):
        if data[i][j][-1] == ",":
            data[i][j] = data[i][j][:-1]
        data[i][j] = fr(data[i][j])
    c = np.array(data[i], dtype = fr).reshape((-1, 1))
    if boo:
        c *= fr("-1")
    c = np.concatenate((c, np.full((cnt, 1), fr("0"), dtype = fr)), axis = 0, dtype = fr)
    input_file.close()
    eye_mat = np.full((m,m), fr("0"), dtype = fr)
    for new_i in range(m):
        eye_mat[new_i, new_i] = fr("1")
    a = np.concatenate((a, eye_mat), axis = 1, dtype = fr)

    answer_dict = {"number_og_cuts" : 0, "initial_solution" : [0]*no_of_vars, "initial_tableau" : None, "final_tableau" : None, "solution_status" : None, "optimal_solution" : [0]*no_of_vars, "optimal_value" : 0}

    c_aux = np.full((n+cnt+m, 1), fr("0"), dtype = fr)
    n_, m_ = n+cnt+m, m 
    tableau = np.full((m_+1, n_+1), fr("0"), dtype = fr)
    init_aux_sol = np.full((n_, 1), fr("0"), dtype = fr)
    for j in range(m_):
        c_aux[n+cnt+j, 0] = fr("1")
        init_aux_sol[n+cnt+j] = b[j]

    tableau[0, 0] = fr("-1") * (c_aux.T @ init_aux_sol)[0,0]
    tableau[0, 1:n_+1] = c_aux.T - np.full((1,m_), fr("1"), dtype = fr)@a
    tableau[1:m_+1, 0] = np.reshape(b, (m,))
    tableau[1:m_+1, 1:n_+1] = a

    def simplex(tableau, basic_vars):
        m, n = tableau.shape
        m -= 1
        n -= 1
        answer = {"solution_status" : "", "optimal_tableau" : None, "basic_variables" : None, "not_none" : 0}
        while True:
            bool_var = 0
            for j in range(1, n+1):
                if (tableau[0, j]) < fr("0"):
                    pivot_col = j
                    break
            else:
                bool_var = 1
            if bool_var == 1:
                break
            theta = None
            pivot_row = None
            for j in range(1, m+1):
                if tableau[j, pivot_col] > fr("0"):
                    if theta == None or tableau[j, 0]/tableau[j, pivot_col] < theta:
                        theta = tableau[j, 0]/tableau[j, pivot_col]
                        pivot_row = j
            if pivot_row == None:
                answer["solution_status"] = "unbounded"
                answer["optimal_tableau"] = tableau
                answer["basic_variables"] = basic_vars
                return answer
            # tableau[pivot_row] /= tableau[pivot_row, pivot_col]
            for j in range(n+1):
                tableau[pivot_row, j] = tableau[pivot_row, j]/tableau[pivot_row, pivot_col]
            for j in range(m+1):
                if j == pivot_row:
                    continue
                for i in range(n+1):
                    tableau[j, i] = tableau[j, i] - tableau[pivot_row, i]*tableau[j, pivot_col]
            basic_vars[pivot_row-1] = pivot_col-1
        answer["solution_status"] = "optimal"
        answer["basic_variables"] = basic_vars
        answer["optimal_tableau"] = tableau
        answer["not_none"] = 1
        return answer

    basic_vars = list(range(n+cnt, n+cnt+m, 1))
    answer_dict["initial_tableau"] = tableau[:, :]
    ans = simplex(tableau, basic_vars)
    if ans["solution_status"] == "unbounded" or (ans["solution_status"] == "optimal" and ans["optimal_tableau"][0, 0] < 0):
        answer_dict["solution_status"] = "infeasible"
    else:
        redundant_constraints = []
        for j in range(m):
            if ans["basic_variables"][j] >= n+cnt:
                for k in range(1, n+cnt+1):
                    if ans["optimal_tableau"][j+1, k] != 0:
                        ans["optimal_tableau"][j+1] /= ans["optimal_tableau"][j+1, k]
                        for z in range(m+1):
                            if z == j+1:
                                continue
                            ans["optimal_tableau"][z] -= ans["optimal_tableau"][j+1]*ans["optimal_tableau"][z, k]
                        ans["basic_variables"][j] = k-1
                        break
                else:
                    redundant_constraints.append(j)
        z = len(redundant_constraints)
        for j in range(z-1, -1, -1):
            ans["basic_variables"].pop(redundant_constraints[j])
            ans["optimal_tableau"] = np.delete(ans["optimal_tableau"], redundant_constraints[j]+1, 0)
        ans["optimal_tableau"] = np.delete(ans["optimal_tableau"], slice(n+cnt,n+cnt+m), axis = 1)

        ans["optimal_tableau"][0, :] *= 0
        for j in range(ans["optimal_tableau"].shape[0] - 1):
            ans["optimal_tableau"][0] += ans["optimal_tableau"][j+1] * c[ans["basic_variables"][j], 0]
        ans["optimal_tableau"][0] = np.concatenate((np.full((1,1), fr("0"), dtype = fr), c.T), axis = 1, dtype = fr) - ans["optimal_tableau"][0]
        
        # answer_dict["initial_tableau"] = ans["optimal_tableau"][:, :]
        ans = simplex(ans["optimal_tableau"], ans["basic_variables"])

        answer_dict["solution_status"] = ans["solution_status"]
        answer_dict["final_tableau"] = ans["optimal_tableau"]
        # answer_dict["final_tableau"] = ans["optimal_tableau"]
        answer_dict["initial_solution"] = [0]*no_of_vars
        if ans["solution_status"] == "optimal":
            for j in range(ans["optimal_tableau"].shape[0] - 1):
                if ans["basic_variables"][j] < n:
                    answer_dict["initial_solution"][ans["basic_variables"][j]] += float(ans["optimal_tableau"][j+1, 0])
    
    ### initial answer created, now moving to integrality
    def dual_simplex(tableau, basic_vars):
        m, n = tableau.shape
        m -= 1
        n -= 1
        answer = {"solution_status" : "", "optimal_tableau" : None, "basic_variables" : None}
        while True:
            bool_var = 0
            for j in range(1, m+1):
                if (tableau[j, 0]) < fr("0"):
                    pivot_row = j
                    break
            else:
                bool_var = 1
            if bool_var == 1:
                break
            theta = None
            pivot_col = None
            for j in range(1, n+1):
                if tableau[pivot_row, j] < fr("0"):
                    if theta == None or tableau[0, j]/abs(tableau[pivot_row, j]) < theta:
                        theta = tableau[0, j]/abs(tableau[pivot_row, j])
                        pivot_col = j
            if pivot_col == None:
                answer["solution_status"] = "infeasible"
                answer["optimal_tableau"] = tableau
                answer["basic_variables"] = basic_vars
                return answer
            for j in range(n+1):
                tableau[pivot_row, j] /= tableau[pivot_row, pivot_col]
            for j in range(m+1):
                if j == pivot_row:
                    continue
                tableau[j] -= tableau[pivot_row]*(tableau[j, pivot_col])
            basic_vars[pivot_row-1] = pivot_col-1
        answer["solution_status"] = "optimal"
        answer["basic_variables"] = basic_vars
        answer["optimal_tableau"] = tableau
        return answer

    if answer_dict["solution_status"] == "optimal":
        while True:
            temp_lst = []
            zzz = len(ans["basic_variables"])
            for i in range(zzz-1, -1, -1):
                if ans["basic_variables"][i] >= no_of_vars+cnt:
                    temp_lst.append(ans["basic_variables"][i])
                    ans["basic_variables"].pop(i)
                    answer_dict["final_tableau"] = np.delete(answer_dict["final_tableau"], i+1, axis = 0)
            temp_lst.sort(reverse=True)
            for i in range(len(temp_lst)):
                answer_dict["final_tableau"] = np.delete(answer_dict["final_tableau"], temp_lst[i]+1, axis = 1)
            # test case 6, 10 left!!
            all_int = 1
            for i in range(answer_dict["final_tableau"].shape[0] - 1):
                if abs(math.floor(answer_dict["final_tableau"][i+1, 0]) - answer_dict["final_tableau"][i+1, 0]) == fr("0"):
                    continue
                else:
                    var = i
                    all_int = 0
                    break
            if all_int == 1:
                break
            else:
                row = np.full((1, answer_dict["final_tableau"].shape[1]), fr("0"), dtype = fr)
                answer_dict["final_tableau"] = np.concatenate((answer_dict["final_tableau"], row), axis=0, dtype = fr)
                for i in range(answer_dict["final_tableau"].shape[1]):
                    answer_dict["final_tableau"][answer_dict["final_tableau"].shape[0] - 1, i] = fr("-1") *answer_dict["final_tableau"][var+1, i] + fr(math.floor(answer_dict["final_tableau"][var+1, i]))
                answer_dict["final_tableau"] = np.concatenate((answer_dict["final_tableau"], np.full((answer_dict["final_tableau"].shape[0], 1), fr("0"), dtype = fr)), axis = 1, dtype = fr)
                answer_dict["final_tableau"][answer_dict["final_tableau"].shape[0] - 1, answer_dict["final_tableau"].shape[1] - 1] = fr("1")
                ans["basic_variables"].append(answer_dict["final_tableau"].shape[1]-1)
                answer_dict["number_og_cuts"] += 1
                ans = dual_simplex(answer_dict["final_tableau"], ans["basic_variables"])
                answer_dict["solution_status"] = ans["solution_status"]
                answer_dict["final_tableau"] = ans["optimal_tableau"]
                if ans["solution_status"] != "optimal":
                    return answer_dict

    answer_dict["optimal_value"] = float(answer_dict["final_tableau"][0, 0]) if boo == 1 else float(fr("-1") * answer_dict["final_tableau"][0, 0])
    answer_dict["optimal_solution"] = [0]*no_of_vars
    ### correct ordering of optimal vector solution
    for j in range(answer_dict["final_tableau"].shape[0] - 1):
        if ans["basic_variables"][j] < no_of_vars:
            answer_dict["optimal_solution"][ans["basic_variables"][j]] += (answer_dict["final_tableau"][j+1, 0])
    return answer_dict

final_ans = gomory_cut_algo()
print("initial_solution: ", end = "")
print(*final_ans["initial_solution"], sep = ", ")
print("final_solution: ", end = "")
print(*final_ans["optimal_solution"], sep = ", ")
print("solution_status:", final_ans["solution_status"])
print("number_of_cuts:", final_ans["number_og_cuts"])
print("optimal_value:", (final_ans["optimal_value"]))