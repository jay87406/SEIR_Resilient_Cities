import base64
import json

import matplotlib.pyplot as pyplot
import networkx
from flask import Flask, jsonify, render_template, request
from flask_cors import CORS
from seirsplus.utilities import *

from models0511 import *
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

    MEAN_INTRACOHORT_DEGREE = int(request_body['MEAN_INTRACOHORT_DEGREE']) or 18
    PCT_CONTACTS_INTERCOHORT = float(request_body['PCT_CONTACTS_INTERCOHORT']) or 0.20
    INIT_EXPOSED = int(request_body['INIT_EXPOSED']) or 30

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
    checkpoints = {'t': [5, 30],
                   'G': [G_quarantine, G_baseline],
                   'p': [0.5 * P_GLOBALINTXN, P_GLOBALINTXN],
                   'theta_E': [0.02, 0.02],
                   'theta_I': [0.02, 0.02],
                   'phi_E': [0.5, 0.5],
                   'phi_I': [0.2, 0.2]}
    T = int(request_body['T']) or 30
    model.run(T=T, checkpoints=checkpoints)
    model.figure_basic(plot_E=True, plot_S=True, plot_R=True, plot_F=False, plot_Q_E=False, plot_Q_I=False, legend=None)
    fig, ax =model.figure_basic(plot_E=True, plot_S=True, plot_R=True, plot_F=False, plot_Q_E=False, plot_Q_I=False, legend=None)
    fig.savefig('test.png')
    img_stream = return_img_stream('test.png')
    return jsonify({
        "img": img_stream
    })

if __name__ == '__main__':
    app.run()
