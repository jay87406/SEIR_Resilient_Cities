from models05112 import *
from networks0217 import *
from sim_loops import *
from seirsplus.utilities import *
import networkx
import matplotlib.pyplot as pyplot
from datetime import datetime

# -------------------------------

# 网络层数，以普陀区为例，有街道、镇、居委、村委设置4
NUM_COHORTS              = 3
# 网络中节点总数(人口) 但设置太多要跑很久，建议等比缩小先看感染曲线趋势
NUM_NODES_PER_COHORT     = 2490
# (街道+镇)、(居委+村委)数量
NUM_TEAMS_PER_COHORT     = [16,215,1078]

MEAN_INTRACOHORT_DEGREE  = 18
PCT_CONTACTS_INTERCOHORT = 0.30

# -------------------------------

N = NUM_NODES_PER_COHORT*NUM_COHORTS

# 模拟开始时当天的感染人数
INIT_EXPOSED = 12

# -------------------------------

G_baseline, cohorts, teams = generate_workplace_contact_network_0218(
                                 num_cohorts=NUM_COHORTS, num_nodes_per_cohort=NUM_NODES_PER_COHORT,
                                 num_teams_per_cohort=NUM_TEAMS_PER_COHORT,
                                 mean_intracohort_degree=MEAN_INTRACOHORT_DEGREE,
                                 pct_contacts_intercohort=PCT_CONTACTS_INTERCOHORT,
                                 farz_params={'alpha':5.0, 'gamma':5.0, 'beta':0.5, 'r':1, 'q':0.0, 'phi':10,
                                              'b':0, 'epsilon':1e-6, 'directed': False, 'weighted': False})

network_info(G_baseline, "Baseline", plot=True)

# --------------------------------

G_quarantine = networkx.classes.function.create_empty_copy(G_baseline)

latentPeriod_mean, latentPeriod_coeffvar = 3.0, 0.6
SIGMA   = 1 / gamma_dist(latentPeriod_mean, latentPeriod_coeffvar, N)
print(SIGMA)

presymptomaticPeriod_mean, presymptomaticPeriod_coeffvar = 2.2, 0.5
LAMDA   = 1 / gamma_dist(presymptomaticPeriod_mean, presymptomaticPeriod_coeffvar, N)
print(LAMDA)

dist_info([1/LAMDA, 1/SIGMA, 1/LAMDA+1/SIGMA], ["latent period", "pre-symptomatic period", "total incubation period"], plot=True, colors=['gold', 'darkorange', 'black'], reverse_plot=True)

symptomaticPeriod_mean, symptomaticPeriod_coeffvar = 4.0, 0.4
GAMMA   = 1 / gamma_dist(symptomaticPeriod_mean, symptomaticPeriod_coeffvar, N)

infectiousPeriod = 1/LAMDA + 1/GAMMA

dist_info([1/LAMDA, 1/GAMMA, 1/LAMDA+1/GAMMA], ["pre-symptomatic period", "(a)symptomatic period", "total infectious period"], plot=True, colors=['darkorange', 'crimson', 'black'], reverse_plot=True)

onsetToHospitalizationPeriod_mean, onsetToHospitalizationPeriod_coeffvar = 11.0, 0.45
ETA     = 1 / gamma_dist(onsetToHospitalizationPeriod_mean, onsetToHospitalizationPeriod_coeffvar, N)

hospitalizationToDischargePeriod_mean, hospitalizationToDischargePeriod_coeffvar = 11.0, 0.45
GAMMA_H = 1 / gamma_dist(hospitalizationToDischargePeriod_mean, hospitalizationToDischargePeriod_coeffvar, N)

dist_info([1/ETA, 1/GAMMA_H, 1/ETA+1/GAMMA_H], ["onset-to-hospitalization period", "hospitalization-to-discharge period", "onset-to-discharge period"], plot=True, colors=['crimson', 'violet', 'black'], reverse_plot=True)

hospitalizationToDeathPeriod_mean, hospitalizationToDeathPeriod_coeffvar = 7.0, 0.45
MU_H    = 1 / gamma_dist(hospitalizationToDeathPeriod_mean, hospitalizationToDeathPeriod_coeffvar, N)

dist_info([1/ETA, 1/MU_H, 1/ETA+1/MU_H], ["onset-to-hospitalization period", "hospitalization-to-death period", "onset-to-death period"], plot=True, colors=['crimson', 'darkgray', 'black'], reverse_plot=True)

# ----------------------------

PCT_ASYMPTOMATIC = 0.25
PCT_HOSPITALIZED = 0.035
PCT_FATALITY = 0.08

# R0_mean 用来设置透过病例数据机算德初的Rt值 (参数名还没敢，叫Rt才是正确的)
# 6.8浦东 3/7    3.35 3/19
R0_mean     = 4.6
R0_coeffvar = 0.2

R0 = gamma_dist(R0_mean, R0_coeffvar, N)

dist_info(R0, "Individual R0", bin_size=0.1, plot=True, colors='crimson')
BETA = 1/infectiousPeriod * R0
print("BETA:",BETA)
P_GLOBALINTXN = 0.2

# ----------------------------

model = SEIRSNetworkModel0511(G=G_baseline, p=P_GLOBALINTXN,
                              beta=BETA, sigma=SIGMA, gamma=GAMMA,
                              G_Q=G_quarantine,
                              initE=INIT_EXPOSED)


# ---------------------------

# 传播模拟结束时间
test_time = '2022/06/25'

# 传播模拟开始时间
start_time = '2022/03/19'
test_DT = datetime.strptime(test_time, '%Y/%m/%d').date()
test_ST = datetime.strptime(start_time, '%Y/%m/%d').date()

# 隔离结束时间
isolation_start_str = '2022/06/30'
isolation_start_date = datetime.strptime(isolation_start_str, '%Y/%m/%d').date()
# 各区开始实施隔离时间
tttttttt_time = '2022/04/01'
start_time_date = datetime.strptime(tttttttt_time, '%Y/%m/%d').date()

isolation_start = int((start_time_date - test_ST).days)
print(isolation_start)
isolation_end = int((isolation_start_date - test_ST).days)
print(isolation_end)
checkpoints = {'t':       [isolation_start, isolation_end],
               'G':       [G_baseline, G_quarantine],
               'p':       [0.9*P_GLOBALINTXN, 0.7*P_GLOBALINTXN],
               'q':       [0.5, 0.5],
               'theta_E': [0.02, 0.02],
               'theta_I': [0.02, 0.02],
               'phi_E':   [0.2, 0.7],
               'phi_I':   [0, 0.15]}

# ---------------------------


T = int((test_DT - test_ST).days)

# --------------------------
# 加入basic模型
SIGMA = SIGMA[0]
GAMMA = GAMMA[0]
MU_I  = MU_H[0]
R0    = R0_mean
BETA  = BETA[1]
model_basic = SEIRSModel(initN   = NUM_NODES_PER_COHORT*10,
                   beta    = BETA,
                   sigma   = SIGMA,
                   gamma   = GAMMA,
                   mu_I    = MU_I,
                   mu_0    = 0,
                   nu      = 0,
                   xi      = 0,
                   beta_Q  = 0.5*(BETA),
                   sigma_Q = SIGMA,
                   gamma_Q = GAMMA,
                   mu_Q    = MU_I,
                   theta_E = 0,
                   theta_I = 0,
                   psi_E   = 1.0,
                   psi_I   = 1.0,
                   initI   = 1,
                   initE   = 1,
                   initQ_E = 0,
                   initQ_I = 0,
                   initR   = 0,
                   initF   = 0)
checkpoints_basic = {'t':       [isolation_start, isolation_end],
               'beta':    [(BETA)*0.3, BETA],
               'theta_E': [0.02, 0.02],
               'theta_I': [0.02, 0.02]
              }
model_basic.run(T=100, checkpoints=checkpoints_basic)
model_basic.figure_infections(vlines=checkpoints_basic['t'])

# -----------------------------

model.run(T, checkpoints=checkpoints)

# -----------------------------

model.figure_basic(plot_S=False, plot_R=True, plot_Q_E=False, plot_Q_I=False, legend=True, vlines=checkpoints['t'], vline_colors=['blue', 'green'], vline_styles=['dashed', 'dotted'], vline_labels=['隔离政策开始', '隔离政策结束'], start_time=start_time, simulation_area='Shanghai')
fig, ax =model.figure_basic(plot_S=False, plot_R=True, plot_Q_E=False, plot_Q_I=False, legend=True, vlines=checkpoints['t'], vline_colors=['blue', 'green'], vline_styles=['dashed', 'dotted'], vline_labels=['隔离政策开始', '隔离政策结束'], start_time=start_time, simulation_area='Shanghai')
fig.savefig('123123.png')