import math
# used for pre-calculate the log values limited by upperLog
maximal_store = 1000000
factorial_num = []

# --------------parameters setting;
# alpha is the significance level;
# cover_rate is the overlapping rate used for remove redundency process;
# alpha_for_edge is the trade-off value between weight and original edge-number, alpha_for_edge is the influence rate of original edge-number to the exact edge-number and thus 1-alpha_for_edge is the influence of weight influence to the exact edge-number
# err_diff is the concussion to avoid local optimum but not global optimum: If new_p-value < original_p-value && original_p-value - new p-value > err diff,  accept the update and continue the search procedure
# max_iter is the loop controller to speed up our algorithm

alpha = 0.01
log_alpha = math.log(alpha)
cover_rate = 0.8
alpha_for_edge = 0.9
err_diff = 0
max_iter = 200
item = 0
name_index = {}
index_name = {}