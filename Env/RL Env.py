from pathlib import Path

import gymnasium as gym
from yby.y_control import Y

from yby.y_nrs import Y_unit
globals().update(Y_unit().unit)

##############################################################################

class Y_env(gym.Env, Y):

    '''
    Paper 2 Env.
    '''

    def __init__(self, verbose=False, **kargs):
        Y.__init__(self)
        for i, v in locals().items():
            setattr(self, i, v)

        self.name = 'Y_env'
        self.start()


    def input_pro(self):

        import numpy as np
        from gymnasium.spaces import Box

        self.cartdim = np.array(self.cartdim)
        self.points = np.array(self.points)
        self.action_dim = self.action_low.__len__()

        # construct action space
        if not self.norm_action:
            self.action_space = Box(low=np.array(self.action_low, dtype=np.float32),
                                    high=np.array(self.action_high, dtype=np.float32),
                                    shape=(self.action_dim,),
                                    dtype=np.float32)
        else:
            self.action_space = Box(low=np.zeros(self.action_dim, dtype=np.float32),
                                    high=np.ones(self.action_dim, dtype=np.float32),
                                    shape=(self.action_dim,),
                                    dtype=np.float32)

        # construct obs space
        obs_low = np.repeat(np.array(self.obs_low), np.prod(np.array(self.cartdim)))
        obs_high = np.repeat(np.array(self.obs_high), np.prod(np.array(self.cartdim)))
        self.observation_space = Box(low=np.array(obs_low, dtype=np.float32),
                                     high=np.array(obs_high, dtype=np.float32),
                                     shape=obs_low.shape,
                                     dtype=np.float32)


    def get_obs(self):

        '''
        Get obs from file.mat
        '''

        import numpy as np
        from scipy.io import loadmat
        from yby.y_mrst import Y_mrst

        root = Path(self.root_tem)

        if self.id_step == 0:
            pres = loadmat(root / 'pres.mat')['pres'][:np.prod(self.cartdim)].flatten()
            sat = loadmat(root / 'sat.mat')['sat'][:np.prod(self.cartdim), 1].flatten()
        else:
            G = root / 'G.mat'

            G = Y_mrst(save=False, verbose=False).mat2dict(G)['G']
            pres = root / 'pres.mat'
            pres = loadmat(pres)['pres'].flatten()
            pres = Y_mrst(save=False, verbose=False).data2df(pres, self.cartdim, frac=True, G=G)
            pres = np.array(pres['data']).flatten()
            sat = root / 'sat.mat'
            sat = loadmat(sat)['sat'][:, 1].flatten()
            sat = Y_mrst(save=False, verbose=False).data2df(sat, self.cartdim, frac=True, G=G)
            sat = np.array(sat['data']).flatten()

        obs = np.append(pres, sat)

        return obs


    def get_volume(self):

        '''
        Get frac volume
        '''

        import numpy as np

        length = self.length[-1]
        height = self.height[-1]
        width = self.width[-1]

        volume = [length, height, width]
        volume = np.prod(volume)

        return volume


    def get_prod(self):

        import numpy as np
        from scipy.io import loadmat
        from yby.y_nrs import Y_unit

        root = Path(self.root_tem)
        prod = loadmat(root / 'qOs.mat')['qOs']
        prod = Y_unit().convert_to(prod, meter ** 3 / day)
        prod = np.abs(prod)

        return prod


    def reset(self, *, seed=2, options=None):

        import numpy as np
        from gymnasium.utils import seeding
        from yby.y_data import Y_krige

        if seed != None:
            self._np_random, seed = seeding.np_random(seed)

        self.id_step = 0
        self.index = []
        self.length = []
        self.height = []
        self.width = []
        self.poro = []
        self.perm = []

        # Y_krige generate file.mat
        bound = {'lb_pres': self.obs_low[0],
                 'ub_pres': self.obs_high[0],
                 'lb_sat': self.obs_low[1],
                 'ub_sat': self.obs_high[1],
                 'sat_ratio': self.sat_ratio}
        cla = Y_krige(self.cartdim, self.points, key=['pres', 'sat'], seed=self.seed, save=False, verbose=False,
                      **bound)
        cla.save_mat(self.root_tem, time=False)

        obs = self.get_obs()
        info = {}
        self.monitor(self.monitor_step, self.id_step, obs)

        return (obs, info)


    def step(self, action):

        import numpy as np
        from yby.y_mrst import Y_mrst
        from yby.y_data import Y_norm

        self.id_step += 1

        # renorm action
        if self.norm_action:
            action = Y_norm(save=False, verbose=False).mm_anti(np.array(action), np.array(self.action_low),
                                                               np.array(self.action_high))

        self.index.append(round(action[0]))
        self.length.append(action[1])
        self.height.append(action[2])
        self.width.append(action[3])
        self.poro.append(action[4])
        self.perm.append(action[5])

        revise = '''
        the_input.frac.index = {};
        the_input.frac.length = {};
        the_input.frac.height = {};
        the_input.frac.width = {};
        the_input.frac.poro = {};
        the_input.frac.perm = {};
        '''.format(self.index, self.length, self.height, self.width, self.poro, self.perm)
        print(revise)

        # calculate obs
        cla = Y_mrst(save=False, verbose=False)
        cla.simu(startup=self.startup, temp_input=self.temp_input, revise=revise, temp_simu=self.temp_simu)
        obs = self.get_obs()

        rwd_action = 0
        volume = self.get_volume()
        rate = self.get_prod()
        prod = rate * self.prod_day
        frac_npv = volume * self.frac_price * self.frac_factor
        prod_npv = prod * self.prod_price * self.prod_factor
        reward = prod_npv - frac_npv + rwd_action

        terminated = False
        truncated = False
        if self.id_step >= self.max_episode_steps:
            terminated = True
        info = {}

        self.monitor(self.monitor_step, self.id_step, obs, index=self.index, action=action, volume=volume, frac_npv=frac_npv,
                     prod=prod, prod_npv=prod_npv, reward_action=rwd_action, reward=reward)

        return (obs, reward, terminated, truncated, info)


    @staticmethod
    def monitor(file, id_step, obs, action=None, index=None, volume=None, frac_npv=None,
                prod=None, prod_npv=None, reward_action=None, reward=None):

        #import sys
        import pickle
        from copy import deepcopy
        #import numpy as np
        import pandas as pd

        #np.set_printoptions(threshold=sys.maxsize)

        obs = deepcopy(obs)
        obs = obs.reshape(2,-1)
        pres = obs[0]
        sat = obs[1]

        data = {'id_step':id_step,
                'pres':pres,
                'sat':sat,
                'action':action,
                'index':index,
                'volume':volume,
                'frac_npv':frac_npv,
                'prod':prod,
                'prod_npv':prod_npv,
                'reward_action':reward_action,
                'reward':reward}

        with open(file, 'ab') as f:
            pickle.dump(data, f)