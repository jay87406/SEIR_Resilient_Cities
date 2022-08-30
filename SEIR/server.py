import base64
import json

import matplotlib.pyplot as pyplot
import networkx
from flask import Flask, jsonify, render_template, request
from flask_cors import CORS
from seirsplus.utilities import *

from models05112 import *
from networks0217 import *
from sim_loops import *

app = Flask(__name__)
# app.debug = True  # Flask内置了调试模式，可以自动重载代码并显示调试信息
CORS(app, supports_credentials=True)

def return_img_stream(img_local_path):
    """
    工具函数:
    获取本地图片流
    :param img_local_path:文件单张图片的本地绝对路径
    :return: 图片流
    """
    import base64
    img_stream = ''
    with open(img_local_path, 'rb') as img_f:
        img_stream = img_f.read()
        img_stream = base64.b64encode(img_stream).decode()
    return img_stream

@app.route("/", methods=["POST"])
def index():
    #############################################
    request_body = json.loads(request.data)
    print(request_body, 99)

    #############################################
    NUM_COHORTS = 2
    NUM_NODES_PER_COHORT = 2490
    # NUM_TEAMS_PER_COHORT     = 16
    NUM_TEAMS_PER_COHORT = [16, 215, 1078]

    N = NUM_NODES_PER_COHORT * NUM_COHORTS

    #-------------------------------输入

    simulation_area = str(request_body['simulation_area']) or 'Shanghai'

    MEAN_INTRACOHORT_DEGREE = int(request_body['MEAN_INTRACOHORT_DEGREE']) or 18
    PCT_CONTACTS_INTERCOHORT = float(request_body['PCT_CONTACTS_INTERCOHORT']) or 0.20
    INIT_EXPOSED = int(request_body['INIT_EXPOSED']) or 30

    # 增加隔离时间 与 开始时间
    isolation_start_str = str(request_body['isolation_start']) or '2022/03/30'
    isolation_end_str = str(request_body['isolation_end']) or '2022/04/30'
    start_time = str(request_body['start_time']) or '2022/03/25'

    # 修改原120行 T = int(request_body['T']) or 30 变成str
    end_time = str(request_body['end_time']) or '2022/04/30'

    #-------------------------------加入日期计算功能
    from datetime import datetime

    isolation_start_date = datetime.strptime(isolation_start_str, '%Y/%m/%d').date()
    isolation_end_date = datetime.strptime(isolation_end_str, '%Y/%m/%d').date()
    start_time_date = datetime.strptime(start_time, '%Y/%m/%d').date()
    end_time_date = datetime.strptime(end_time, '%Y/%m/%d').date()

    # 隔离开始时间计算
    isolation_start = int((isolation_start_date - start_time_date).days)
    # 隔离结束时间计算
    isolation_end = int((isolation_end_date - start_time_date).days)
    # 模拟时间 T 计算
    T = int((end_time_date - start_time_date).days)

    #-------------------------------
    G_baseline, cohorts, teams = generate_workplace_contact_network_0218(
        num_cohorts=NUM_COHORTS, num_nodes_per_cohort=NUM_NODES_PER_COHORT,
        num_teams_per_cohort=NUM_TEAMS_PER_COHORT,
        mean_intracohort_degree=MEAN_INTRACOHORT_DEGREE,
        pct_contacts_intercohort=PCT_CONTACTS_INTERCOHORT,
        farz_params={'alpha': 5.0, 'gamma': 5.0, 'beta': 0.5, 'r': 1, 'q': 0.0, 'phi': 10,
                     'b': 0, 'epsilon': 1e-6, 'directed': False, 'weighted': False})
    G_quarantine = networkx.classes.function.create_empty_copy(G_baseline)
    #--------------------------------------------
    latentPeriod_mean, latentPeriod_coeffvar = 3.0, 0.6
    SIGMA = 1 / gamma_dist(latentPeriod_mean, latentPeriod_coeffvar, N)
    presymptomaticPeriod_mean, presymptomaticPeriod_coeffvar = 2.2, 0.5
    LAMDA = 1 / gamma_dist(presymptomaticPeriod_mean, presymptomaticPeriod_coeffvar, N)
    symptomaticPeriod_mean, symptomaticPeriod_coeffvar = 4.0, 0.4
    GAMMA = 1 / gamma_dist(symptomaticPeriod_mean, symptomaticPeriod_coeffvar, N)
    infectiousPeriod = 1 / LAMDA + 1 / GAMMA
    onsetToHospitalizationPeriod_mean, onsetToHospitalizationPeriod_coeffvar = 11.0, 0.45
    ETA = 1 / gamma_dist(onsetToHospitalizationPeriod_mean, onsetToHospitalizationPeriod_coeffvar, N)
    hospitalizationToDischargePeriod_mean, hospitalizationToDischargePeriod_coeffvar = 11.0, 0.45
    GAMMA_H = 1 / gamma_dist(hospitalizationToDischargePeriod_mean, hospitalizationToDischargePeriod_coeffvar, N)
    hospitalizationToDeathPeriod_mean, hospitalizationToDeathPeriod_coeffvar = 7.0, 0.45
    MU_H = 1 / gamma_dist(hospitalizationToDeathPeriod_mean, hospitalizationToDeathPeriod_coeffvar, N)

    #--------------------------------------------
    PCT_ASYMPTOMATIC = 0.25
    PCT_HOSPITALIZED = 0.035
    PCT_FATALITY = 0.08

    R0 = float(request_body['R0']) or 2.5

    BETA = 1 / infectiousPeriod * R0
    P_GLOBALINTXN = 0.2
    #############################################
    model = SEIRSNetworkModel0511(G=G_baseline, p=P_GLOBALINTXN,
                              beta=BETA, sigma=SIGMA, gamma=GAMMA,
                              G_Q=G_quarantine,
                              initE=INIT_EXPOSED)
    checkpoints = {'t': [isolation_start, isolation_end],
                   'G': [G_quarantine, G_baseline],
                   'p': [0.5 * P_GLOBALINTXN, P_GLOBALINTXN],
                   'theta_E': [0.02, 0.02],
                   'theta_I': [0.02, 0.02],
                   'phi_E': [0.5, 0.5],
                   'phi_I': [0.2, 0.2]}
    #T = int(request_body['T']) or 30

    # --------------------------
    # 加入basic模型
    SIGMA = SIGMA[0]
    GAMMA = GAMMA[0]
    MU_I = MU_H[0]
    R0 = R0
    BETA = BETA[0]
    model_basic = SEIRSModel(initN=NUM_NODES_PER_COHORT * 10,
                             beta=BETA,
                             sigma=SIGMA,
                             gamma=GAMMA,
                             mu_I=MU_I,
                             mu_0=0,
                             nu=0,
                             xi=0,
                             beta_Q=0.5 * (BETA),
                             sigma_Q=SIGMA,
                             gamma_Q=GAMMA,
                             mu_Q=MU_I,
                             theta_E=0,
                             theta_I=0,
                             psi_E=1.0,
                             psi_I=1.0,
                             initI=1,
                             initE=0,
                             initQ_E=0,
                             initQ_I=0,
                             initR=0,
                             initF=0)
    checkpoints_basic = {'t': [isolation_start, isolation_end],
                         'beta': [(BETA) * 0.3, BETA],
                         'theta_E': [0.02, 0.02],
                         'theta_I': [0.02, 0.02]
                         }
    model_basic.run(T=100, checkpoints=checkpoints_basic)
    model_basic.figure_infections(vlines=checkpoints_basic['t'])

    # -----------------------------

    model.run(T=T, checkpoints=checkpoints)

    model.figure_basic(plot_S=False, plot_R=True, plot_Q_E=False, plot_Q_I=False, legend=True, vlines=checkpoints['t'],
                       vline_colors=['blue', 'green'], vline_styles=['dashed', 'dotted'],
                       vline_labels=['隔离政策开始', '隔离政策结束'], start_time=start_time, simulation_area='Shanghai')
    fig, ax = model.figure_basic(plot_S=False, plot_R=True, plot_Q_E=False, plot_Q_I=False, legend=True,
                                 vlines=checkpoints['t'], vline_colors=['blue', 'green'],
                                 vline_styles=['dashed', 'dotted'], vline_labels=['隔离政策开始', '隔离政策结束'],
                                 start_time=start_time, simulation_area='Shanghai')
    fig.savefig('test.png')
    img_stream = return_img_stream('test.png')
    return jsonify({
        "img": img_stream
    })

if __name__ == '__main__':
    app.run()
