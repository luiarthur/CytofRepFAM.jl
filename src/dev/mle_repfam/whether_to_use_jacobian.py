import torch
from torch import nn
from torch.distributions import Exponential, Gamma
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

# NOTE:
# I wasn't sure whether I should be using the Jacobian in MAP estimation when I
# have continuous parameters that don't have support on the entire real line.
# Because I'm optimizing, and not sampling, I don't need a jacobian. See the
# results by running this script.
#
# Sketch of proof that a Jacobian isn't needed:
# Say we are interested in solving for theta* = argmin_theta f(theta), with
# theta > 0.  This is the same as solving for gamma* = argmin_gamma g(gamma),
# where g(gamma) = f(exp(gamma)), for unconstrained gamma. If gamma* minimizes
# g, then exp(gamma*) minimizes f.  That is, theta* = exp(gamma*). No Jacobian
# needed here.

class Model(nn.Module):
    def __init__(self, y, prior_shape=2, prior_rate=3, use_jacobian=False):
        super(Model, self).__init__()

        self.y = y
        self.y_sum = y.sum()
        self.n = y.shape[0]
        self.prior_shape = prior_shape
        self.prior_rate = prior_rate
        self.use_jacobian = use_jacobian
        
        # self.log_theta = torch.nn.Parameter(torch.zeros([]))
        self.log_theta = torch.nn.Parameter(torch.ones([]))
        self.loss_hist = []

    def logprob(self, log_theta):
        theta = log_theta.exp()
        if self.use_jacobian:
            log_jacobian = log_theta
        else:
            log_jacobian = 0

        return ((n + self.prior_shape - 1) * log_theta - 
                theta * (self.y_sum + self.prior_rate) +
                log_jacobian)

    def forward(self):
        return self.logprob(self.log_theta)
    
    def fit(self, niter=10000, tol=1e-6, lr=.0001, verbose=1):
        opt = torch.optim.Adam(self.parameters(), lr=lr)

        for t in range(niter):
            # Update model
            loss = -self.forward() / self.n
            opt.zero_grad()
            loss.backward()
            opt.step()
            self.loss_hist.append(loss.item())

            # iterator.set_description('Loss: {:.2f}'.format(loss.item()))

            if (t > 10 and
                abs(self.loss_hist[-1] / self.loss_hist[-2] - 1) < tol and
                verbose):
                print('Convergence detected!')
                break



if __name__ == '__main__':
    torch.manual_seed(0)
    theta_true = 5.3
    # NOTE: see n=3, n=10, n=30, n=1000
    # For small n, the MAP estimate with jacobian is far from truth.  As n
    # increases, the difference is small, because the posterior gets more
    # peaked.
    n = 5
    y = Exponential(theta_true).sample((n, ))

    model = Model(y) 
    model.fit(lr=.01)

    model_with_jacobian = Model(y, use_jacobian=True)
    model_with_jacobian.fit(lr=.01)

    true_map_est = ((n + model.prior_shape - 1) / 
                    (model.prior_rate + model.y_sum))

    map_est = model.log_theta.exp()
    map_est_with_jacobian = model_with_jacobian.log_theta.exp()

    print('True MAP est: {}'.format(true_map_est))
    print('MAP est (w/o jacobian): {}'.format(map_est))
    print('MAP est (w/ jacobian): {}'.format(map_est_with_jacobian))

    # NOTE: The correct way to do MAP estimation is to NOT use a jacobian!

    plt.plot(model_with_jacobian.loss_hist, label='w/ jacobian')
    plt.plot(model.loss_hist, label='w/o jacobian')
    plt.legend()
    plt.show()

    # Plot results
    post_shape = model.n + model.prior_shape
    post_rate = model.prior_rate + model.y_sum
    # plt.hist(Gamma(post_shape, post_rate).sample((100000, )), bins=64)
    x = torch.arange(0.0, 10.0, step=.01)
    d = model.logprob(x.log())
    plt.plot(x, d, lw=2)
    plt.axvline(true_map_est, color='green', lw=2,
                label='MAP estimate (true)', ls='--')
    plt.axvline(map_est, color='orange', lw=2, label='MAP estimate (no jacobian)',
                ls=':')
    plt.axvline(map_est_with_jacobian, color='blue', lw=2,
                label='MAP estimate (with jacobian)')
    plt.xlabel('theta')
    plt.ylabel('log density')
    plt.legend()
    plt.show()
