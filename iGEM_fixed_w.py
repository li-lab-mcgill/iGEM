import torch
import scipy.io
import numpy as np
from tqdm import tqdm, trange
import pickle
import pandas as pd
from joblib import parallel_backend
        
class iGEM_fixed_w(object):
    def __init__(self, X1,X2,aa_mat,ab,bb, params, fixed_model_w):
        #dup for convenience
        self.params = params
        params_dict = params
        self.X1 = X1
        self.X2 = X2
        
        aa = aa_mat
        
        self.k = params_dict['k']
        self.max_iter = params_dict['max_iter']
        self.tol = params_dict['tol']
        pi_aa = params_dict['pi_aa']
        pi_ab = params_dict['pi_ab']
        pi_bb = params_dict['pi_bb']
        pi_xr = params_dict['pi_xr']
        pi_xc = params_dict['pi_xc']
        self.xc_alpha1_coef = params_dict['xc_alpha1_coef']
        self.xc_h2_coef = params_dict['xc_h2_coef']
        
        self.pat = X1.shape[0]
        
        #use pi for threshold
        self.pi_aa = pi_aa
        self.pi_ab = pi_ab
        self.pi_bb = pi_bb
        self.pi_xr = pi_xr
        self.pi_xc = pi_xc
        self.aa = aa
        self.ab = ab
        
        self.omic_num = params_dict['omic_num']
            
        torch.manual_seed(2020)
            
        self.W = fixed_model_w
        self.loss = []
        
        self.use_alpha = params_dict['use_alpha']
        self.use_poisson = params_dict['use_poisson']
        
        
        if self.omic_num == 2:
            if self.use_alpha[0]:
                self.rho1 = params_dict['rho1']
                self.alpha1 = torch.ones(self.k, self.rho1.shape[0])
                self.H1 = torch.mm(self.alpha1, self.rho1)
            else:
                self.H1 = torch.ones(self.k, X1.shape[1])
            if self.use_alpha[1]:
                self.rho2 = params_dict['rho2']
                self.alpha2 = torch.ones(self.k, self.rho2.shape[0])
                self.H2 = torch.mm(self.alpha2, self.rho2)
            else:
                self.H2 = torch.ones(self.k, X2.shape[1])
        else:
            if self.use_alpha[0]:
                self.rho1 = params_dict['rho1']
                self.alpha1 = torch.ones(self.k, self.rho1.shape[0])
                self.H1 = torch.mm(self.alpha1, self.rho1)
            else:
                self.H1 = torch.ones(self.k, X1.shape[1])
     
        if aa==None:
            self.aa = torch.zeros(X1.shape[1],X1.shape[1])
            
#line 91     
    def train(self):

        aa = self.aa
        ab = self.ab
        denom_min = 1e-12
               
        with parallel_backend('threading', n_jobs=-1):
            with torch.no_grad():
                prev_loss = 0
                
                for iters in trange(self.max_iter):
                    if self.omic_num == 1:
                        #update W
                        #self.W = self.W * (torch.mm(self.X1, self.H1.t())) / (torch.mm(self.W, (torch.mm(self.H1, self.H1.t()) + self.pi_xr*torch.eye(self.k)))+denom_min)
                        #update H
                        if not self.use_poisson[0]:
                            if self.use_alpha[0]:
                                self.alpha1 = self.alpha1 * torch.mm(torch.mm(self.W.t(), self.X1),self.rho1.t()) / (torch.mm(torch.mm((torch.mm(self.W.t(), self.W)), self.H1),self.rho1.t()) + self.pi_xc*torch.mm(torch.eye(self.k, self.k),self.alpha1) + denom_min)
                                self.H1 = self.H1 * (torch.mm(self.W.t(), self.X1)) / (torch.mm((torch.mm(self.W.t(), self.W)), self.H1) + self.pi_xc*torch.mm(torch.eye(self.k, self.k),self.H1) + denom_min)
                            else:
                                self.H1 = self.H1 * (torch.mm(self.W.t(), self.X1) + self.pi_aa*torch.mm(self.H1, aa)) / (torch.mm((torch.mm(self.W.t(), self.W) + self.pi_xc*torch.ones(self.k)), self.H1) + denom_min)
                        else:
                            if self.use_alpha[0]:
                                self.alpha1 = self.alpha1 * torch.mm((torch.mm(self.W.t(), (self.X1 / torch.mm(self.W, self.H1))) + self.pi_aa*torch.mm(self.H1, aa) ), self.rho1.t())/((torch.mm(torch.mm(self.W.t(), torch.ones(self.X1.shape[0], self.X1.shape[1])), self.rho1.t()) + self.pi_xc*torch.mm(torch.eye(self.k, self.k),self.alpha1)) + denom_min)
                                self.H1 = torch.mm(self.alpha1, self.rho1)
                            else:
                                self.H1 = self.H1 * (torch.mm(self.W.t(), (self.X1 / torch.mm(self.W, self.H1))) + self.pi_aa*torch.mm(self.H1, aa) )/ (torch.mm(self.W.t(), torch.ones(self.X1.shape[0], self.X1.shape[1])) + self.pi_xc*torch.mm(torch.eye(self.k, self.k),self.H1) + denom_min)
                    elif self.omic_num == 2:
                        #self.W = self.W * (torch.mm(self.X1, self.H1.t()) + torch.mm(self.X2, self.H2.t())) / (torch.mm(self.W, (torch.mm(self.H1, self.H1.t()) + torch.mm(self.H2, self.H2.t()) + self.pi_xr*torch.eye(self.k)))+denom_min)
                        if not self.use_poisson[1]:
                            if self.use_alpha[0]:
                                tmp_H1 = self.H1 * (torch.mm(self.W.t(), self.X1) + torch.mm(self.alpha1, self.rho1) + self.pi_aa*torch.mm(self.H1, aa) + 0.5 * self.pi_ab*torch.mm(self.H2, ab.t()))/(torch.mm((torch.mm(self.W.t(), self.W) + self.pi_xc*torch.ones(self.k, self.k)), self.H1) + self.H1 + denom_min)
                                self.alpha1 = self.alpha1 * torch.mm(self.H1, self.rho1.t()) / (torch.mm((torch.mm(self.alpha1, self.rho1)), self.rho1.t()) + self.xc_alpha1_coef*torch.mm(torch.ones(self.k, self.k),self.alpha1) + denom_min)

                            else:
                                tmp_H1 = self.H1 * (torch.mm(self.W.t(), self.X1) + self.pi_aa*torch.mm(self.H1, aa) + 0.5*self.pi_ab*torch.mm(self.H2, ab.t())) / (torch.mm((torch.mm(self.W.t(), self.W) + self.pi_xc*torch.ones(self.k)), self.H1) + denom_min)

                            if self.use_alpha[1]:
                                self.H2 = torch.mm(self.alpha2, self.rho2)
                            else:
                                self.H2 = self.H2 * (torch.mm(self.W.t(), self.X2) + 0.5*self.pi_ab*torch.mm(self.H1, ab)) / (torch.mm((torch.mm(self.W.t(), self.W) + self.xc_h2_coef*self.pi_xc*torch.ones(self.k)), self.H2) + denom_min)
                            self.H1 = tmp_H1
                        
                        else:
                            if self.use_alpha[0]:
                                self.alpha1 = self.alpha1 * torch.mm((torch.mm(self.W.t(), (self.X1 / torch.mm(self.W, self.H1))) + self.pi_aa*torch.mm(self.H1, aa) ), self.rho1.t())/((torch.mm(torch.mm(self.W.t(), torch.ones(self.X1.shape[0], self.X1.shape[1])), self.rho1.t()) + self.pi_xc*torch.mm(torch.eye(self.k, self.k),self.alpha1)) + 1e-8)
                                self.H1 = torch.mm(self.alpha1, self.rho1)
                            else:
                                self.H1 = self.H1 * (torch.mm(self.W.t(), self.X1) + self.pi_aa*torch.mm(self.H1, aa) + 0.5*self.pi_ab*torch.mm(self.H2, ab.t())) / (torch.mm((torch.mm(self.W.t(), self.W) + self.pi_xc*torch.ones(self.k)), self.H1) + 1e-8)

                            if self.use_alpha[1]:
                                self.H2 = torch.mm(self.alpha2, self.rho2)
                            else:
                                self.H2 = self.H2 * (torch.mm(self.W.t(), self.X2)) / (torch.mm((torch.mm(self.W.t(), self.W)), self.H2) + 1e-8)
                            #self.H1 = tmp_H1
                    
                    #self.H1 = tmp_H1
                    
                    if self.omic_num == 2:
                        if torch.isnan(self.H1).any():
                            print(f'H1 fails at {iters}')
                            break
                        if torch.isnan(self.H2).any():
                            print(f'H2 fails at {iters}')
                            break
                    else:
                        if torch.isnan(self.H1).any():
                            print(f'H1 fails at {iters}')
                            break
                    
                    if self.use_alpha[0]:
                        self.prev_terms = [self.alpha1, self.H1, self.W]
                    else:
                        self.prev_terms = [self.H1, self.W]
        
                    #SE
                    if self.omic_num == 1:
                        err1 = self.X1 - torch.mm(self.W, self.H1)
                        err1 = torch.pow(err1, 2).sum()
                        err = err1
                        err2 = 0
                        err3 = 0
                    elif self.omic_num == 2:
                        err1 = self.X1 - torch.mm(self.W, self.H1)
                        err1 = torch.pow(err1, 2).sum()
                        if self.use_alpha[0]:
                            err1 += torch.pow((self.H1 - torch.mm(self.alpha1, self.rho1)), 2).sum()
                        err2 = self.X2 - torch.mm(self.W, self.H2)
                        err2 = torch.pow(err2, 2).sum()
                        err = err1 + err2
                        err3 = 0                    

                    trr = 0
                    if self.pi_aa != 0:
                        trr += -self.pi_aa*torch.trace(torch.mm(torch.mm(self.H1, aa), self.H1.t()))
                    if self.pi_ab != 0 and self.omic_num > 1:
                        trr += -self.pi_ab*torch.trace(torch.mm(torch.mm(self.H1, ab), self.H2.t()))
                    reg_xr = self.pi_xr*(self.W*self.W).sum()
                    xc = 0
                    if self.omic_num == 2:
                        if self.use_alpha[0]:
                            xc_alpha = self.xc_alpha1_coef * (self.alpha1*self.alpha1).sum()
                            xc_h1 = self.pi_xc * (self.H1*self.H1).sum()
                            xc_tmp = xc_alpha + xc_h1
                        else:
                            xc_tmp = (self.H1*self.H1).sum()
                        xc += xc_tmp
                        if self.use_alpha[1]:
                            xc_tmp = (self.alpha2*self.alpha2).sum()
                            xc_tmp += (self.H2*self.H2).sum()
                        else:
                            xc_h2 = self.xc_h2_coef * (self.H2*self.H2).sum()
                            xc_tmp = xc_h2
                        xc += xc_tmp
                    else:
                        if self.use_alpha[0]:
                            xc_tmp = (self.alpha1*self.alpha1).sum()
                        else:
                            xc_tmp = (self.H1*self.H1).sum()
                        xc += xc_tmp
                       

                    reg = reg_xr + xc

                    new_loss = err + trr + reg
                    loss_lst = []
                    
                    if iters % 100 == 0:
                        print(f'Loss = {new_loss}')

                    if abs(new_loss - prev_loss) < self.tol:
                        self.loss.append(loss_lst)
                        print('converged')
                        break
                    else:
                        self.loss.append(loss_lst)
                        prev_loss = new_loss