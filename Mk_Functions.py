
import numpy as np
import math

def zone111(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    k = np.ones_like(beta) * (np.sqrt((rho_t + rho_c)**2 * n**2 + (((-2 * xi + 2) * alpha + 2 * xi) * rho_c + 2 * rho_t * (1 + (xi - 1) * alpha)) * n + xi) - 1 + (-rho_c - rho_t) * n) / (xi - 1)
    M = -3 * E * b * beta * epsilon_cr * ((xi / 3 - 1 / 3) * k**3 + (1 + (rho_t + rho_c) * n) * k**2 + (-1 + ((2 * alpha - 2) * rho_c - 2 * rho_t * alpha) * n) * k + 1 / 3 + ((alpha - 1)**2 * rho_c + rho_t * alpha**2) * n) * h**2 / (3 * k - 3) 
    return k, M




def zone211(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    # Calculate `k`
    k = np.real((np.sqrt(
        ((-2 * n * alpha * (-(-1 + beta)**2 * eta_1 + beta**2 * xi - 2 * beta + 1) * (rho_c - rho_t) + 
                    2 * (n * rho_t + xi / 2) * (-1 + beta)**2 * eta_1 + 
                    n * (n * rho_c**2 + (2 * n * rho_t + 2 * xi) * rho_c + n * rho_t**2) * beta**2 + 
                    (4 * n * rho_t + 2 * xi) * beta - 2 * n * rho_t - xi) * beta**2)) + 
         (-eta_1 + (-rho_c - rho_t) * n) * beta**2 + (2 * eta_1 - 2) * beta - eta_1 + 1)/ ((xi - eta_1) * beta**2 + (2 * eta_1 - 2) * beta - eta_1 + 1))

    # Calculate `M`
    M = -epsilon_cr * E * b * h**2 * (
        (((-eta_1 / 3 + xi / 3) * k**3 + (eta_1 + (rho_c + rho_t) * n) * k**2 + 
          (-eta_1 + ((2 * alpha - 2) * rho_c - 2 * rho_t * alpha) * n) * k + eta_1 / 3 + 
          ((alpha**2 - 2 * alpha + 1) * rho_c + rho_t * alpha**2) * n) * beta**3 + 
         (k - 1)**3 * (eta_1 - 1) * beta**2 / 2 - (k - 1)**3 * (eta_1 - 1) / 6) / 
        (beta**2 * (k - 1))
    )
    
    return k, M


def zone212(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    k = np.real((-(-1 + beta)**2 * eta_1 - n * (eta_s * rho_t + rho_c) * beta**2 + 
                (kappa * n * eta_s * rho_t - kappa * n * rho_t - 2) * beta + 
                np.sqrt(beta**2 * (((eta_s * rho_t + rho_c) * beta - rho_t * kappa * (eta_s - 1))**2 * n**2 + 
                ((2 * ((-alpha + 1) * eta_1 + xi * alpha) * eta_s * rho_t - 
                  2 * rho_c * (-alpha * eta_1 + xi * (alpha - 1))) * beta**2 + 
                ((((4 * alpha - 4) * eta_1 - 2 * xi * kappa - 4 * alpha + 4) * eta_s + 
                  2 * xi * kappa) * rho_t - 4 * alpha * rho_c * (eta_1 - 1)) * beta - 
                2 * ((alpha - 1) * eta_s * rho_t - rho_c * alpha) * (eta_1 - 1)) * n + 
                xi * (beta**2 * eta_1 + (-2 * eta_1 + 2) * beta + eta_1 - 1))) + 1) / 
                (-(-1 + beta)**2 * eta_1 + beta**2 * xi - 2 * beta + 1))

    M = (h**2 * E * epsilon_cr * (((eta_1 - xi) * k**3 + ((-3 * eta_s * rho_t - 3 * rho_c) * n - 3 * eta_1) * k**2 + 
         ((6 * eta_s * rho_t * alpha - 6 * rho_c * (alpha - 1)) * n + 3 * eta_1) * k + 
         (-3 * eta_s * rho_t * alpha**2 - 3 * rho_c * (alpha - 1)**2) * n - eta_1) * beta**3 - 
         3 * ((eta_1 - 1) * k**2 + (-2 * rho_t * kappa * (eta_s - 1) * n - 2 * eta_1 + 2) * k + 
         2 * rho_t * alpha * kappa * (eta_s - 1) * n + eta_1 - 1) * (k - 1) * beta**2 / 2 + 
         (k - 1)**3 * (eta_1 - 1) / 2) * b / (3 * beta**2 * (k - 1)))

    return k, M



def zone221(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    k = np.real(((-(-1 + beta)**2 * eta_1 + xi * omega * (beta + omega) * eta_c - n * (rho_c + rho_t) * beta**2 + 
                 (-omega * xi - 2) * beta - xi * omega**2 + 
                 np.sqrt(beta**2 * (beta**2 * (rho_c + rho_t)**2 * n**2 + 
                 ((-2 * eta_c * (-rho_t * alpha + rho_c * (alpha - 1)) * xi + 
                   2 * eta_1 * ((-alpha + 1) * rho_t + rho_c * alpha)) * beta**2 + 
                 (-4 * (alpha - 0.5) * omega * (eta_c - 1) * (rho_c - rho_t) * xi - 
                   4 * (eta_1 - 1) * ((-alpha + 1) * rho_t + rho_c * alpha)) * beta - 
                 2 * ((-alpha + 1) * rho_t + rho_c * alpha) * (omega**2 * (eta_c - 1) * xi - eta_1 + 1)) * n - 
                 xi * (-eta_1 * eta_c * beta**2 + 2 * eta_c * (eta_1 - 1) * beta + 
                       omega**2 * (eta_c - 1) * xi - eta_c * (eta_1 - 1)))) + 1) / 
                 (xi * (beta + omega)**2 * eta_c - (-1 + beta)**2 * eta_1 - 
                  2 * beta * omega * xi - xi * omega**2 - 2 * beta + 1)))

    M = (b * ((((-xi * eta_c + eta_1) * beta**3 + (-(3 * omega * (eta_c - 1) * xi) / 2 + 
             3 / 2 - (3 * eta_1) / 2) * beta**2 + omega**3 * (eta_c - 1) * xi / 2 - 
             1 / 2 + eta_1 / 2) * k**3 + ((-3 * eta_1 + (-3 * rho_c - 3 * rho_t) * n) * 
             beta**3 + ((3 * omega * (eta_c - 1) * xi) / 2 - 9 / 2 + (9 * eta_1) / 2) * 
             beta**2 - (3 * omega**3 * (eta_c - 1) * xi) / 2 + 3 / 2 - (3 * eta_1) / 2) * 
             k**2 + ((3 * eta_1 + ((-6 * alpha + 6) * rho_c + 6 * rho_t * alpha) * n) * 
             beta**3 + (9 / 2 - (9 * eta_1) / 2) * beta**2 + (3 * omega**3 * (eta_c - 1) * 
             xi) / 2 - 3 / 2 + (3 * eta_1) / 2) * k + (-eta_1 + ((-3 * alpha**2 + 6 * alpha - 
             3) * rho_c - 3 * rho_t * alpha**2) * n) * beta**3 + (-3 / 2 + (3 * eta_1) / 2) * 
             beta**2 - omega**3 * (eta_c - 1) * xi / 2 + 1 / 2 - eta_1 / 2) * 
             E * epsilon_cr * h**2 / (3 * (k - 1) * beta**2)))
    
    return k, M




def zone222(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    k = np.real((xi * omega * (beta + omega) * eta_c - (-1 + beta)**2 * eta_1 - xi * omega**2 - beta * omega * xi - 
                n * (eta_s * rho_t + rho_c) * beta**2 + (kappa * n * eta_s * rho_t - kappa * n * rho_t - 2) * beta + 
                np.sqrt((((beta * eta_s - kappa * (eta_s - 1)) * rho_t + rho_c * beta)**2 * n**2 + 
                ((2 * (xi * alpha * eta_c - (alpha - 1) * eta_1) * eta_s * beta**2 + 
                ((((4 * alpha * omega - 2 * kappa - 2 * omega) * eta_c - 4 * alpha * omega + 2 * omega) * eta_s + 
                2 * kappa * eta_c) * xi + 4 * eta_s * (eta_1 - 1) * (alpha - 1)) * beta + 
                2 * omega * ((alpha * omega - kappa - omega) * eta_s + kappa) * (eta_c - 1) * xi - 
                2 * eta_s * (eta_1 - 1) * (alpha - 1)) * rho_t - 
                2 * ((eta_c * (alpha - 1) * xi - eta_1 * alpha) * beta**2 + 
                (2 * omega * (alpha - 0.5) * (eta_c - 1) * xi + 2 * alpha * (eta_1 - 1)) * beta + 
                alpha * (omega**2 * (eta_c - 1) * xi - eta_1 + 1)) * rho_c) * n - 
                xi * (-eta_1 * eta_c * beta**2 + 2 * eta_c * (eta_1 - 1) * beta + omega**2 * (eta_c - 1) * xi - eta_c * (eta_1 - 1))) * beta**2) + 1) / 
                (xi * (beta + omega)**2 * eta_c - xi * omega**2 - 2 * beta * omega * xi - (-1 + beta)**2 * eta_1 - 2 * beta + 1))

    M = ((((-2 * xi * eta_c + 2 * eta_1) * beta**3 + (-3 * omega * (eta_c - 1) * xi - 3 * eta_1 + 3) * beta**2 + 
         omega**3 * (eta_c - 1) * xi + eta_1 - 1) * k**3 + 
         (((-6 * eta_s * rho_t - 6 * rho_c) * n - 6 * eta_1) * beta**3 + 
         (6 * rho_t * kappa * (eta_s - 1) * n + 3 * omega * (eta_c - 1) * xi + 9 * eta_1 - 9) * beta**2 - 
         3 * omega**3 * (eta_c - 1) * xi - 3 * eta_1 + 3) * k**2 + 
         (((12 * eta_s * rho_t * alpha - 12 * rho_c * (alpha - 1)) * n + 6 * eta_1) * beta**3 + 
         (-6 * rho_t * kappa * (alpha + 1) * (eta_s - 1) * n - 9 * eta_1 + 9) * beta**2 + 
         3 * omega**3 * (eta_c - 1) * xi + 3 * eta_1 - 3) * k + 
         ((-6 * eta_s * rho_t * alpha**2 - 6 * rho_c * (alpha - 1)**2) * n - 2 * eta_1) * beta**3 + 
         (6 * rho_t * alpha * kappa * (eta_s - 1) * n + 3 * eta_1 - 3) * beta**2 - omega**3 * (eta_c - 1) * xi - eta_1 + 1) * 
         b * E * epsilon_cr * h**2 / (6 * (k - 1) * beta**2))

    return k, M




def zone311(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    k = np.real((np.sqrt(beta**2 * (4 * (n * (rho_c - rho_t) * alpha + n * rho_t + xi / 2) * 
                (-1 + beta_1) * (-beta_1 / 2 + beta - 1 / 2) * eta_1 + 
                2 * (beta - beta_1)**2 * (n * (rho_c - rho_t) * alpha + n * rho_t + xi / 2) * eta_2 - 
                2 * n * (rho_c - rho_t) * (beta**2 * xi - 2 * beta + 1) * alpha + 
                (n * rho_c**2 + (2 * n * rho_t + 2 * xi) * rho_c + n * rho_t**2) * n * beta**2 + 
                (4 * n * rho_t + 2 * xi) * beta - 2 * n * rho_t - xi)) + 
                (-eta_2 + (-rho_c - rho_t) * n) * beta**2 + 
                ((-2 * eta_1 + 2 * eta_2) * beta_1 + 2 * eta_1 - 2) * beta + 
                (-eta_2 + eta_1) * beta_1**2 - eta_1 + 1) / 
                ((xi - eta_2) * beta**2 + ((-2 * eta_1 + 2 * eta_2) * beta_1 + 2 * eta_1 - 2) * beta + 
                 (-eta_2 + eta_1) * beta_1**2 - eta_1 + 1))

    M = (h**2 * b * (((2 * eta_2 - 2 * xi) * beta**3 + ((-3 * eta_2 + 3 * eta_1) * beta_1 - 3 * eta_1 + 3) * beta**2 + 
         (eta_2 - eta_1) * beta_1**3 + eta_1 - 1) * k**3 + 
         ((-6 * eta_2 + (-6 * rho_c - 6 * rho_t) * n) * beta**3 + ((9 * eta_2 - 9 * eta_1) * beta_1 + 9 * eta_1 - 9) * beta**2 + 
         (-3 * eta_2 + 3 * eta_1) * beta_1**3 - 3 * eta_1 + 3) * k**2 + 
         ((6 * eta_2 + ((-12 * alpha + 12) * rho_c + 12 * rho_t * alpha) * n) * beta**3 + ((-9 * eta_2 + 9 * eta_1) * beta_1 - 9 * eta_1 + 9) * 
         beta**2 + (3 * eta_2 - 3 * eta_1) * beta_1**3 + 3 * eta_1 - 3) * k + 
         (-2 * eta_2 + ((-6 * alpha**2 + 12 * alpha - 6) * rho_c - 6 * rho_t * alpha**2) * n) * beta**3 + 
         ((3 * eta_2 - 3 * eta_1) * beta_1 + 3 * eta_1 - 3) * beta**2 + (-eta_2 + eta_1) * beta_1**3 - eta_1 + 1) * 
         epsilon_cr * E / (6 * beta**2 * (k - 1)))

    return k, M




def zone312(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    k = np.real(((-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1 - (beta_1 - beta)**2 * eta_2 - 
                 n * (eta_s * rho_t + rho_c) * beta**2 + (kappa * n * eta_s * rho_t - kappa * n * rho_t - 2) * beta + 
                 np.sqrt(beta**2 * (((eta_s * rho_t + rho_c) * beta - rho_t * kappa * (eta_s - 1))**2 * n**2 + 
                 ((2 * ((xi - eta_2) * alpha + eta_2) * eta_s * rho_t - 2 * rho_c * ((xi - eta_2) * alpha - xi)) * beta**2 + 
                 (((((-4 * eta_1 + 4 * eta_2) * beta_1 + 4 * eta_1 - 4) * alpha + (4 * eta_1 - 4 * eta_2) * beta_1 - 
                 2 * xi * kappa - 4 * eta_1 + 4) * eta_s + 2 * xi * kappa) * rho_t + 
                 4 * alpha * rho_c * ((eta_1 - eta_2) * beta_1 - eta_1 + 1)) * beta + 
                 2 * ((eta_1 - eta_2) * beta_1**2 - eta_1 + 1) * (eta_s * (alpha - 1) * rho_t - rho_c * alpha)) * n - 
                 xi * (-beta**2 * eta_2 + ((-2 * eta_1 + 2 * eta_2) * beta_1 + 2 * eta_1 - 2) * beta + 
                       (eta_1 - eta_2) * beta_1**2 - eta_1 + 1))) + 1) / 
                 ((-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1 - (beta_1 - beta)**2 * eta_2 + beta**2 * xi - 2 * beta + 1))

    M = (-b * E * epsilon_cr * h**2 * (((-2 * eta_2 + 2 * xi) * beta**3 + ((-3 * eta_1 + 3 * eta_2) * beta_1 + 3 * eta_1 - 3) * beta**2 + 
          (eta_1 - eta_2) * beta_1**3 - eta_1 + 1) * k**3 + 
          (((6 * eta_s * rho_t + 6 * rho_c) * n + 6 * eta_2) * beta**3 + (-6 * rho_t * kappa * (eta_s - 1) * n + 
          (9 * eta_1 - 9 * eta_2) * beta_1 - 9 * eta_1 + 9) * beta**2 + (-3 * eta_1 + 3 * eta_2) * beta_1**3 + 3 * eta_1 - 3) * k**2 + 
          (((-12 * eta_s * alpha * rho_t + 12 * rho_c * (alpha - 1)) * n - 6 * eta_2) * beta**3 + 
          (6 * rho_t * kappa * (alpha + 1) * (eta_s - 1) * n + (-9 * eta_1 + 9 * eta_2) * beta_1 + 9 * eta_1 - 9) * beta**2 + 
          (3 * eta_1 - 3 * eta_2) * beta_1**3 - 3 * eta_1 + 3) * k + 
          ((6 * eta_s * rho_t * alpha**2 + 6 * rho_c * (alpha - 1)**2) * n + 2 * eta_2) * beta**3 + 
          (-6 * rho_t * alpha * kappa * (eta_s - 1) * n + (3 * eta_1 - 3 * eta_2) * beta_1 - 3 * eta_1 + 3) * beta**2 + 
          (-eta_1 + eta_2) * beta_1**3 + eta_1 - 1) / (6 * beta**2 * (k - 1)))

    return k, M




def zone321(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    k = np.real((np.sqrt((2 * (-1 + beta_1) * (beta_1 - 2 * beta + 1) * (-xi * eta_c / 2 + ((rho_t - rho_c) * alpha - rho_t) * n) * eta_1 + 
                 2 * xi * ((beta_1 - beta)**2 * eta_2 / 2 - xi * omega**2 / 2 + n * (beta + omega)**2 * (rho_t - rho_c) * alpha - 
                 n * rho_t * omega**2 - n * beta * (rho_t - rho_c) * omega + rho_c * n * beta**2 + beta - 1 / 2) * eta_c - 
                 2 * (beta_1 - beta)**2 * ((rho_t - rho_c) * alpha - rho_t) * n * eta_2 + xi**2 * omega**2 - 
                 4 * omega * ((rho_t - rho_c) * (beta + omega / 2) * alpha - rho_t * omega / 2 - beta * (rho_t - rho_c) / 2) * n * xi + 
                 (-4 * (rho_t - rho_c) * (beta - 1 / 2) * alpha + n * (rho_c + rho_t)**2 * beta**2 + 4 * beta * rho_t - 2 * rho_t) * n) * beta**2) + 
                 (-eta_2 + (-rho_t - rho_c) * n) * beta**2 + ((eta_c - 1) * omega * xi + (-2 * eta_1 + 2 * eta_2) * beta_1 + 2 * eta_1 - 2) * beta + 
                 omega**2 * (eta_c - 1) * xi + (eta_1 - eta_2) * beta_1**2 - eta_1 + 1) / 
                 ((xi * eta_c - eta_2) * beta**2 + (2 * (eta_c - 1) * omega * xi + (-2 * eta_1 + 2 * eta_2) * beta_1 + 2 * eta_1 - 2) * beta + 
                  omega**2 * (eta_c - 1) * xi + (eta_1 - eta_2) * beta_1**2 - eta_1 + 1))

    M = -(((2 * xi * eta_c - 2 * eta_2) * beta**3 + ((-3 * eta_1 + 3 * eta_2) * beta_1 + 3 * eta_1 - 3 + (3 * eta_c - 3) * omega * xi) * beta**2 + 
          (eta_1 - eta_2) * beta_1**3 - eta_1 + 1 + (-eta_c + 1) * omega**3 * xi) * k**3 + 
          ((6 * eta_2 + (6 * rho_c + 6 * rho_t) * n) * beta**3 + ((9 * eta_1 - 9 * eta_2) * beta_1 - 9 * eta_1 + 9 + (-3 * eta_c + 3) * omega * xi) * 
          beta**2 + (-3 * eta_1 + 3 * eta_2) * beta_1**3 + 3 * eta_1 - 3 + (3 * eta_c - 3) * omega**3 * xi) * k**2 + 
          ((-6 * eta_2 + ((12 * alpha - 12) * rho_c - 12 * rho_t * alpha) * n) * beta**3 + ((-9 * eta_1 + 9 * eta_2) * beta_1 + 9 * eta_1 - 9) * 
          beta**2 + (3 * eta_1 - 3 * eta_2) * beta_1**3 - 3 * eta_1 + 3 + (-3 * eta_c + 3) * omega**3 * xi) * k + 
          (2 * eta_2 + ((6 * alpha**2 - 12 * alpha + 6) * rho_c + 6 * rho_t * alpha**2) * n) * beta**3 + 
          ((3 * eta_1 - 3 * eta_2) * beta_1 - 3 * eta_1 + 3) * beta**2 + (-eta_1 + eta_2) * beta_1**3 + eta_1 - 1 + 
          (eta_c - 1) * omega**3 * xi) * h**2 * b * epsilon_cr * E / (6 * (k - 1) * beta**2)

    return k, M



def zone322(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    k = np.real(((-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1 + xi * omega * (beta + omega) * eta_c - 
                 (beta_1 - beta)**2 * eta_2 - xi * omega * (beta + omega) - n * (eta_s * rho_t + rho_c) * beta**2 + 
                 (kappa * n * eta_s * rho_t - kappa * n * rho_t - 2) * beta + 
                 np.sqrt(beta**2 * (((eta_s * beta - kappa * (eta_s - 1)) * rho_t + beta * rho_c)**2 * n**2 + 
                 ((2 * (xi * alpha * eta_c - eta_2 * (alpha - 1)) * eta_s * beta**2 + 
                 (((4 * omega * (eta_c - 1) * alpha + (-2 * kappa - 2 * omega) * eta_c + 2 * omega) * xi - 
                 4 * ((eta_1 - eta_2) * beta_1 - eta_1 + 1) * (alpha - 1)) * eta_s + 2 * xi * kappa * eta_c) * beta + 
                 (2 * omega * (eta_c - 1) * (alpha * omega - kappa - omega) * xi + 
                 2 * ((eta_1 - eta_2) * beta_1**2 - eta_1 + 1) * (alpha - 1)) * eta_s + 
                 2 * xi * kappa * omega * (eta_c - 1)) * rho_t - 
                 2 * ((eta_c * (alpha - 1) * xi - alpha * eta_2) * beta**2 + 
                 (2 * (eta_c - 1) * (alpha - 0.5) * omega * xi - 
                 2 * ((eta_1 - eta_2) * beta_1 - eta_1 + 1) * alpha) * beta + 
                 alpha * (omega**2 * (eta_c - 1) * xi + (eta_1 - eta_2) * beta_1**2 - eta_1 + 1)) * rho_c) * n - 
                 (-eta_2 * eta_c * beta**2 - 2 * eta_c * ((eta_1 - eta_2) * beta_1 - eta_1 + 1) * beta + 
                 omega**2 * (eta_c - 1) * xi + eta_c * ((eta_1 - eta_2) * beta_1**2 - eta_1 + 1)) * xi)) + 1) / 
                 ((-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1 + xi * (beta + omega)**2 * eta_c - 
                 (beta_1 - beta)**2 * eta_2 + (-2 * beta * omega - omega**2) * xi - 2 * beta + 1))

    M = (epsilon_cr * h**2 * (((-2 * xi * eta_c + 2 * eta_2) * beta**3 + 
          ((3 * eta_1 - 3 * eta_2) * beta_1 - 3 * eta_1 + 3 + (-3 * eta_c + 3) * omega * xi) * beta**2 + 
          (-eta_1 + eta_2) * beta_1**3 + eta_1 - 1 + (eta_c - 1) * omega**3 * xi) * k**3 + 
          (((-6 * eta_s * rho_t - 6 * rho_c) * n - 6 * eta_2) * beta**3 + 
          (6 * rho_t * kappa * (eta_s - 1) * n + (-9 * eta_1 + 9 * eta_2) * beta_1 + 9 * eta_1 - 9 + 
          (3 * eta_c - 3) * omega * xi) * beta**2 + 
          (3 * eta_1 - 3 * eta_2) * beta_1**3 - 3 * eta_1 + 3 + (-3 * eta_c + 3) * omega**3 * xi) * k**2 + 
          (((12 * eta_s * rho_t * alpha - 12 * rho_c * (alpha - 1)) * n + 6 * eta_2) * beta**3 + 
          (-6 * rho_t * kappa * (alpha + 1) * (eta_s - 1) * n + (9 * eta_1 - 9 * eta_2) * beta_1 - 9 * eta_1 + 9) * beta**2 + 
          (-3 * eta_1 + 3 * eta_2) * beta_1**3 + 3 * eta_1 - 3 + (3 * eta_c - 3) * omega**3 * xi) * k + 
          ((-6 * eta_s * rho_t * alpha**2 - 6 * rho_c * (alpha - 1)**2) * n - 2 * eta_2) * beta**3 + 
          (6 * rho_t * alpha * kappa * (eta_s - 1) * n + (-3 * eta_1 + 3 * eta_2) * beta_1 + 3 * eta_1 - 3) * beta**2 + 
          (eta_1 - eta_2) * beta_1**3 - eta_1 + 1 + (-eta_c + 1) * omega**3 * xi) * E * b / (6 * (k - 1) * beta**2))

    return k, M



def zone411(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    k = np.real((np.sqrt(beta**2 * (-2 * (n * (rho_c - rho_t) * alpha + n * rho_t + xi / 2) * (-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1 + 
                 2 * (n * (rho_c - rho_t) * alpha + n * rho_t + xi / 2) * (beta_1 - beta_2) * (beta_1 - 2 * beta + beta_2) * eta_2 + 
                 2 * (n * (rho_c - rho_t) * alpha + n * rho_t + xi / 2) * (beta_2 - beta)**2 * eta_3 - 
                 2 * n * (rho_c - rho_t) * (beta**2 * xi - 2 * beta + 1) * alpha + 
                 (n * rho_c**2 + (2 * n * rho_t + 2 * xi) * rho_c + n * rho_t**2) * n * beta**2 + 
                 (4 * n * rho_t + 2 * xi) * beta - 2 * n * rho_t - xi)) + 
                 (-eta_3 + (-rho_c - rho_t) * n) * beta**2 + 
                 ((-2 * eta_1 + 2 * eta_2) * beta_1 + (-2 * eta_2 + 2 * eta_3) * beta_2 + 2 * eta_1 - 2) * beta + 
                 (-eta_2 + eta_1) * beta_1**2 + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1) / 
                 ((xi - eta_3) * beta**2 + ((-2 * eta_1 + 2 * eta_2) * beta_1 + (-2 * eta_2 + 2 * eta_3) * beta_2 + 
                 2 * eta_1 - 2) * beta + (-eta_2 + eta_1) * beta_1**2 + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1))

    M = (epsilon_cr * h**2 * (((2 * eta_3 - 2 * xi) * beta**3 + ((-3 * eta_2 + 3 * eta_1) * beta_1 + (3 * eta_2 - 3 * eta_3) * beta_2 - 
          3 * eta_1 + 3) * beta**2 + (eta_2 - eta_1) * beta_1**3 + (-eta_2 + eta_3) * beta_2**3 + eta_1 - 1) * k**3 + 
          ((-6 * eta_3 + (-6 * rho_c - 6 * rho_t) * n) * beta**3 + ((9 * eta_2 - 9 * eta_1) * beta_1 + (-9 * eta_2 + 9 * eta_3) * 
          beta_2 + 9 * eta_1 - 9) * beta**2 + (-3 * eta_2 + 3 * eta_1) * beta_1**3 + (3 * eta_2 - 3 * eta_3) * beta_2**3 - 
          3 * eta_1 + 3) * k**2 + ((6 * eta_3 + ((-12 * alpha + 12) * rho_c + 12 * rho_t * alpha) * n) * beta**3 + 
          ((-9 * eta_2 + 9 * eta_1) * beta_1 + (9 * eta_2 - 9 * eta_3) * beta_2 - 9 * eta_1 + 9) * beta**2 + 
          (3 * eta_2 - 3 * eta_1) * beta_1**3 + (-3 * eta_2 + 3 * eta_3) * beta_2**3 + 3 * eta_1 - 3) * k + 
          (-2 * eta_3 + ((-6 * alpha**2 + 12 * alpha - 6) * rho_c - 6 * rho_t * alpha**2) * n) * beta**3 + 
          ((3 * eta_2 - 3 * eta_1) * beta_1 + (-3 * eta_2 + 3 * eta_3) * beta_2 + 3 * eta_1 - 3) * beta**2 + 
          (-eta_2 + eta_1) * beta_1**3 + (eta_2 - eta_3) * beta_2**3 - eta_1 + 1) * b * E / (6 * beta**2 * (k - 1)))

    return k, M


'''

def zone412(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    k = np.real(((-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1 - (beta_1 - beta_2) * (beta_1 - 2 * beta + beta_2) * eta_2 - 
                 (beta_2 - beta)**2 * eta_3 - n * (eta_s * rho_t + rho_c) * beta**2 + 
                 (kappa * n * eta_s * rho_t - kappa * n * rho_t - 2) * beta + 
                 np.sqrt(beta**2 * (((eta_s * rho_t + rho_c) * beta - rho_t * kappa * (eta_s - 1))**2 * n**2 + 
                 ((2 * eta_s * ((xi - eta_3) * alpha + eta_3) * rho_t - 2 * rho_c * ((xi - eta_3) * alpha - xi)) * beta**2 + 
                 (((((-4 * eta_1 + 4 * eta_2) * beta_1 + (-4 * eta_2 + 4 * eta_3) * beta_2 + 4 * eta_1 - 4) * alpha - 
                 2 * xi * kappa + (4 * eta_1 - 4 * eta_2) * beta_1 + (4 * eta_2 - 4 * eta_3) * beta_2 - 
                 4 * eta_1 + 4) * eta_s + 2 * xi * kappa) * rho_t + 
                 4 * ((eta_1 - eta_2) * beta_1 + (eta_2 - eta_3) * beta_2 - eta_1 + 1) * alpha * rho_c) * beta - 
                 2 * (-eta_s * (alpha - 1) * rho_t + rho_c * alpha) * ((eta_1 - eta_2) * beta_1**2 + 
                 (eta_2 - eta_3) * beta_2**2 - eta_1 + 1)) * n - 
                 (-beta**2 * eta_3 + ((-2 * eta_1 + 2 * eta_2) * beta_1 + (-2 * eta_2 + 2 * eta_3) * beta_2 + 
                 2 * eta_1 - 2) * beta + (eta_1 - eta_2) * beta_1**2 + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1) * xi)) + 1) / 
                 ((-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1 - (beta_1 - beta_2) * (beta_1 - 2 * beta + beta_2) * eta_2 - 
                 (beta_2 - beta)**2 * eta_3 + beta**2 * xi - 2 * beta + 1))

    M = (-h**2 * epsilon_cr * E * b * (((-2 * eta_3 + 2 * xi) * beta**3 + 
          ((-3 * eta_1 + 3 * eta_2) * beta_1 + (-3 * eta_2 + 3 * eta_3) * beta_2 + 3 * eta_1 - 3) * beta**2 + 
          (eta_1 - eta_2) * beta_1**3 + (eta_2 - eta_3) * beta_2**3 - eta_1 + 1) * k**3 + 
          (((6 * eta_s * rho_t + 6 * rho_c) * n + 6 * eta_3) * beta**3 + 
          (-6 * rho_t * kappa * (eta_s - 1) * n + (9 * eta_1 - 9 * eta_2) * beta_1 + 
          (9 * eta_2 - 9 * eta_3) * beta_2 - 9 * eta_1 + 9) * beta**2 + 
          (-3 * eta_1 + 3 * eta_2) * beta_1**3 + (-3 * eta_2 + 3 * eta_3) * beta_2**3 + 3 * eta_1 - 3) * k**2 + 
          (((-12 * eta_s * rho_t * alpha + 12 * rho_c * (alpha - 1)) * n - 6 * eta_3) * beta**3 + 
          (6 * rho_t * kappa * (alpha + 1) * (eta_s - 1) * n + (-9 * eta_1 + 9 * eta_2) * beta_1 + 
          (-9 * eta_2 + 9 * eta_3) * beta_2 + 9 * eta_1 - 9) * beta**2 + 
          (3 * eta_1 - 3 * eta_2) * beta_1**3 + (3 * eta_2 - 3 * eta_3) * beta_2**3 - 3 * eta_1 + 3) * k + 
          ((6 * eta_s * rho_t * alpha**2 + 6 * rho_c * (alpha - 1)**2) * n + 2 * eta_3) * beta**3 + 
          (-6 * rho_t * alpha * kappa * (eta_s - 1) * n + (3 * eta_1 - 3 * eta_2) * beta_1 + 
          (3 * eta_2 - 3 * eta_3) * beta_2 - 3 * eta_1 + 3) * beta**2 + 
          (-eta_1 + eta_2) * beta_1**3 + (-eta_2 + eta_3) * beta_2**3 + eta_1 - 1) / (6 * beta**2 * (k - 1)))

    return k, M
'''
def zone412(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3,
            eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s,
            rho_c, rho_t):
    
    # Compute k
    #
    # First, build the inner expression under the square root.
    sqrt_inner = (
        beta**2 * (
            (
                ((eta_s * rho_t + rho_c) * beta - rho_t * kappa * (eta_s - 1))**2 * n**2
                + (
                    (2 * eta_s * (((xi - eta_3) * alpha) + eta_3) * rho_t
                     - 2 * rho_c * (((xi - eta_3) * alpha) - xi)
                    ) * beta**2
                    + (
                        (
                          (((-4 * eta_1 + 4 * eta_2) * beta_1
                             + (-4 * eta_2 + 4 * eta_3) * beta_2
                             + 4 * eta_1 - 4) * alpha
                           - 2 * xi * kappa
                           + (4 * eta_1 - 4 * eta_2) * beta_1
                           + (4 * eta_2 - 4 * eta_3) * beta_2
                           - 4 * eta_1 + 4)
                        * eta_s
                        + 2 * xi * kappa
                        ) * rho_t
                    + 4 * ((eta_1 - eta_2) * beta_1 + (eta_2 - eta_3) * beta_2 - eta_1 + 1)
                      * alpha * rho_c
                    ) * beta
                - 2 * (-eta_s * (alpha - 1) * rho_t + rho_c * alpha)
                  * ((eta_1 - eta_2) * beta_1**2 + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1)
            ) * n
            - (
                - beta**2 * eta_3
                + (((-2 * eta_1 + 2 * eta_2) * beta_1
                    + (-2 * eta_2 + 2 * eta_3) * beta_2
                    + 2 * eta_1 - 2) * beta)
                + (eta_1 - eta_2) * beta_1**2
                + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1
            ) * xi
        )
    ))
    
    # Now compute the numerator and denominator for k.
    num_k = (
          (-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1
        - (beta_1 - beta_2) * (beta_1 - 2 * beta + beta_2) * eta_2
        - (beta_2 - beta)**2 * eta_3
        - n * (eta_s * rho_t + rho_c) * beta**2
        + (kappa * n * eta_s * rho_t - kappa * n * rho_t - 2) * beta
        + np.sqrt(sqrt_inner)
        + 1
    )
    den_k = (
          (-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1
        - (beta_1 - beta_2) * (beta_1 - 2 * beta + beta_2) * eta_2
        - (beta_2 - beta)**2 * eta_3
        + beta**2 * xi
        - 2 * beta + 1
    )
    
    # Use np.real to remove any (tiny) imaginary parts.
    k = np.real(num_k) / den_k

    # Compute M.
    #
    # It is helpful to break the long expression into four groups.
    term1 = (
          ((-2 * eta_3 + 2 * xi) * beta**3
           + (((-3 * eta_1 + 3 * eta_2) * beta_1
               + (-3 * eta_2 + 3 * eta_3) * beta_2 + 3 * eta_1 - 3)
              * beta**2
              + (eta_1 - eta_2) * beta_1**3
              + (eta_2 - eta_3) * beta_2**3 - eta_1 + 1)
          ) * k**3
    )
    term2 = (
          ((((6 * eta_s * rho_t + 6 * rho_c) * n + 6 * eta_3) * beta**3)
           + ((-6 * rho_t * kappa * (eta_s - 1) * n
               + (9 * eta_1 - 9 * eta_2) * beta_1
               + (9 * eta_2 - 9 * eta_3) * beta_2 - 9 * eta_1 + 9)
              * beta**2)
           + (-3 * eta_1 + 3 * eta_2) * beta_1**3
           + (-3 * eta_2 + 3 * eta_3) * beta_2**3 + 3 * eta_1 - 3
          ) * k**2
    )
    term3 = (
          (((( -12 * eta_s * rho_t * alpha + 12 * rho_c * (alpha - 1)) * n - 6 * eta_3) * beta**3)
           + ((6 * rho_t * kappa * (alpha + 1) * (eta_s - 1) * n
               + (-9 * eta_1 + 9 * eta_2) * beta_1
               + (-9 * eta_2 + 9 * eta_3) * beta_2 + 9 * eta_1 - 9)
              * beta**2)
           + (3 * eta_1 - 3 * eta_2) * beta_1**3
           + (3 * eta_2 - 3 * eta_3) * beta_2**3 - 3 * eta_1 + 3
          ) * k
    )
    term4 = (
          ((6 * eta_s * rho_t * alpha**2 + 6 * rho_c * (alpha - 1)**2) * n + 2 * eta_3) * beta**3
          + (-6 * rho_t * alpha * kappa * (eta_s - 1) * n
             + (3 * eta_1 - 3 * eta_2) * beta_1
             + (3 * eta_2 - 3 * eta_3) * beta_2 - 3 * eta_1 + 3) * beta**2
          + (-eta_1 + eta_2) * beta_1**3
          + (-eta_2 + eta_3) * beta_2**3 + eta_1 - 1
    )
    
    M = (-h**2 * epsilon_cr * E * b *
         (term1 + term2 + term3 + term4)
         / (6 * beta**2 * (k - 1))
    )
    
    return k, M   

def zone421(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    k = np.real((np.sqrt((-2 * (-1 + beta_1) * (beta_1 - 2 * beta + 1) * (xi * eta_c / 2 + ((rho_c - rho_t) * alpha + rho_t) * n) * eta_1 + 
                 2 * (xi * eta_c / 2 + ((rho_c - rho_t) * alpha + rho_t) * n) * (beta_1 - 2 * beta + beta_2) * (beta_1 - beta_2) * eta_2 - 
                 2 * xi * (-(beta_2 - beta)**2 * eta_3 / 2 + xi * omega**2 / 2 + n * (beta + omega)**2 * (rho_c - rho_t) * alpha + 
                 n * rho_t * omega**2 - n * beta * (rho_c - rho_t) * omega - rho_c * n * beta**2 - beta + 1 / 2) * eta_c + 
                 2 * ((rho_c - rho_t) * alpha + rho_t) * n * (beta_2 - beta)**2 * eta_3 + omega**2 * xi**2 + 
                 4 * n * omega * ((beta + omega / 2) * (rho_c - rho_t) * alpha + rho_t * omega / 2 - beta * (rho_c - rho_t) / 2) * xi + 
                 n * (4 * (beta - 1 / 2) * (rho_c - rho_t) * alpha + n * (rho_c + rho_t)**2 * beta**2 + 4 * beta * rho_t - 2 * rho_t)) * beta**2) + 
                 (-eta_3 + (-rho_c - rho_t) * n) * beta**2 + 
                 (omega * (eta_c - 1) * xi + (-2 * eta_1 + 2 * eta_2) * beta_1 + (-2 * eta_2 + 2 * eta_3) * beta_2 + 2 * eta_1 - 2) * beta + 
                 omega**2 * (eta_c - 1) * xi + (eta_1 - eta_2) * beta_1**2 + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1) / 
                 ((xi * eta_c - eta_3) * beta**2 + (2 * omega * (eta_c - 1) * xi + (-2 * eta_1 + 2 * eta_2) * beta_1 + 
                 (-2 * eta_2 + 2 * eta_3) * beta_2 + 2 * eta_1 - 2) * beta + omega**2 * (eta_c - 1) * xi + 
                 (eta_1 - eta_2) * beta_1**2 + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1))

    M = (-E * epsilon_cr * (((2 * xi * eta_c - 2 * eta_3) * beta**3 + 
          ((-3 * eta_1 + 3 * eta_2) * beta_1 + (-3 * eta_2 + 3 * eta_3) * beta_2 + 3 * eta_1 - 3 + (3 * eta_c - 3) * omega * xi) * beta**2 + 
          (eta_1 - eta_2) * beta_1**3 + (eta_2 - eta_3) * beta_2**3 - eta_1 + 1 + (-eta_c + 1) * omega**3 * xi) * k**3 + 
          ((6 * eta_3 + (6 * rho_c + 6 * rho_t) * n) * beta**3 + 
          ((9 * eta_1 - 9 * eta_2) * beta_1 + (9 * eta_2 - 9 * eta_3) * beta_2 - 9 * eta_1 + 9 + (-3 * eta_c + 3) * omega * xi) * beta**2 + 
          (-3 * eta_1 + 3 * eta_2) * beta_1**3 + (-3 * eta_2 + 3 * eta_3) * beta_2**3 + 3 * eta_1 - 3 + (3 * eta_c - 3) * omega**3 * xi) * k**2 + 
          ((-6 * eta_3 + ((12 * alpha - 12) * rho_c - 12 * rho_t * alpha) * n) * beta**3 + 
          ((-9 * eta_1 + 9 * eta_2) * beta_1 + (-9 * eta_2 + 9 * eta_3) * beta_2 + 9 * eta_1 - 9) * beta**2 + 
          (3 * eta_1 - 3 * eta_2) * beta_1**3 + (3 * eta_2 - 3 * eta_3) * beta_2**3 - 3 * eta_1 + 3 + (-3 * eta_c + 3) * omega**3 * xi) * k + 
          (2 * eta_3 + ((6 * alpha**2 - 12 * alpha + 6) * rho_c + 6 * rho_t * alpha**2) * n) * beta**3 + 
          ((3 * eta_1 - 3 * eta_2) * beta_1 + (3 * eta_2 - 3 * eta_3) * beta_2 - 3 * eta_1 + 3) * beta**2 + 
          (-eta_1 + eta_2) * beta_1**3 + (-eta_2 + eta_3) * beta_2**3 + eta_1 - 1 + (eta_c - 1) * omega**3 * xi) * b * h**2 / (6 * (k - 1) * beta**2))

    return k, M

def zone422(beta_z4, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    beta = beta_z4  # Ensure beta_z4 is a NumPy array for vectorization
    
    beta_z4 = np.array(beta_z4)
    k = np.zeros(len(beta_z4))
    for i in range(len(beta_z4)):
        beta = beta_z4[i]  # Changed from beta_z4[i,0] to beta_z4[i]
        k[i] = ((-1 + beta_1)*(beta_1 - 2*beta + 1)*eta_1 - (beta_1 - beta_2)*(beta_1 - 2*beta + beta_2)*eta_2 + 
                xi*omega*(beta + omega)*eta_c - (beta_2 - beta)**2*eta_3 - xi*omega*(beta + omega) - 
                n*(eta_s*rho_t + rho_c)*beta**2 + (kappa*n*eta_s*rho_t - kappa*n*rho_t - 2)*beta + 
                np.sqrt((((eta_s*beta - kappa*(eta_s - 1))*rho_t + beta*rho_c)**2*n**2 + 
                        ((2*eta_s*(eta_c*xi*alpha - eta_3*(alpha - 1))*beta**2 + 
                          (((4*omega*(eta_c - 1)*alpha + (-2*kappa - 2*omega)*eta_c + 2*omega)*xi - 
                            4*(alpha - 1)*((eta_1 - eta_2)*beta_1 + (eta_2 - eta_3)*beta_2 - eta_1 + 1))*eta_s + 
                           2*eta_c*xi*kappa)*beta + 
                          (2*omega*(alpha*omega - kappa - omega)*(eta_c - 1)*xi + 
                           2*(alpha - 1)*((eta_1 - eta_2)*beta_1**2 + (eta_2 - eta_3)*beta_2**2 - eta_1 + 1))*eta_s + 
                          2*xi*kappa*omega*(eta_c - 1))*rho_t - 
                         2*rho_c*((eta_c*(alpha - 1)*xi - alpha*eta_3)*beta**2 + 
                                 (2*(alpha - 0.5)*omega*(eta_c - 1)*xi - 
                                  2*((eta_1 - eta_2)*beta_1 + (eta_2 - eta_3)*beta_2 - eta_1 + 1)*alpha)*beta + 
                                 alpha*(omega**2*(eta_c - 1)*xi + (eta_1 - eta_2)*beta_1**2 + 
                                       (eta_2 - eta_3)*beta_2**2 - eta_1 + 1)))*n - 
                        xi*(-eta_c*eta_3*beta**2 - 
                            2*eta_c*((eta_1 - eta_2)*beta_1 + (eta_2 - eta_3)*beta_2 - eta_1 + 1)*beta + 
                            omega**2*(eta_c - 1)*xi + 
                            ((eta_1 - eta_2)*beta_1**2 + (eta_2 - eta_3)*beta_2**2 - eta_1 + 1)*eta_c))*beta**2) + 
                1)/((-1 + beta_1)*(beta_1 - 2*beta + 1)*eta_1 - 
                    (beta_1 - beta_2)*(beta_1 - 2*beta + beta_2)*eta_2 + 
                    xi*(beta + omega)**2*eta_c - (beta_2 - beta)**2*eta_3 + 
                    (-2*beta*omega - omega**2)*xi - 2*beta + 1)
    
    beta = beta_z4
    M = E * b * (((-2 * xi * eta_c + 2 * eta_3) * beta**3 + 
                  ((3 * eta_1 - 3 * eta_2) * beta_1 + (3 * eta_2 - 3 * eta_3) * beta_2 - 3 * eta_1 + 3 + 
                   (-3 * eta_c + 3) * omega * xi) * beta**2 + 
                  (-eta_1 + eta_2) * beta_1**3 + (-eta_2 + eta_3) * beta_2**3 + eta_1 - 1 + 
                  (eta_c - 1) * omega**3 * xi) * k**3 + 
                 (((-6 * eta_s * rho_t - 6 * rho_c) * n - 6 * eta_3) * beta**3 + 
                  (6 * rho_t * kappa * (eta_s - 1) * n + (-9 * eta_1 + 9 * eta_2) * beta_1 + 
                   (-9 * eta_2 + 9 * eta_3) * beta_2 + 9 * eta_1 - 9 + 
                   (3 * eta_c - 3) * omega * xi) * beta**2 + 
                  (3 * eta_1 - 3 * eta_2) * beta_1**3 + (3 * eta_2 - 3 * eta_3) * beta_2**3 - 3 * eta_1 + 3 + 
                  (-3 * eta_c + 3) * omega**3 * xi) * k**2 + 
                 (((12 * eta_s * rho_t * alpha - 12 * rho_c * (alpha - 1)) * n + 6 * eta_3) * beta**3 + 
                  (-6 * rho_t * kappa * (alpha + 1) * (eta_s - 1) * n + (9 * eta_1 - 9 * eta_2) * beta_1 + 
                   (9 * eta_2 - 9 * eta_3) * beta_2 - 9 * eta_1 + 9) * beta**2 + 
                  (-3 * eta_1 + 3 * eta_2) * beta_1**3 + (-3 * eta_2 + 3 * eta_3) * beta_2**3 + 3 * eta_1 - 3 + 
                  (3 * eta_c - 3) * omega**3 * xi) * k + 
                 ((-6 * eta_s * rho_t * alpha**2 - 6 * rho_c * (alpha - 1)**2) * n - 2 * eta_3) * beta**3 + 
                 (6 * rho_t * alpha * kappa * (eta_s - 1) * n + (-3 * eta_1 + 3 * eta_2) * beta_1 + 
                  (-3 * eta_2 + 3 * eta_3) * beta_2 + 3 * eta_1 - 3) * beta**2 + 
                 (eta_1 - eta_2) * beta_1**3 + (eta_2 - eta_3) * beta_2**3 - eta_1 + 1 + 
                 (-eta_c + 1) * omega**3 * xi) * h**2 * epsilon_cr / (6 * (k - 1) * beta**2)
    
    return k, M

    
    ''' 
    # Common terms used multiple times
    common_term = (eta_1 - eta_2) * beta_1 + (eta_2 - eta_3) * beta_2 - eta_1 + 1
    beta_squared_terms = (eta_1 - eta_2) * beta_1**2 + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1
    
    numerator_k = (
        ((-1 + beta_1) * (beta_1 - 2*beta + 1) * eta_1) 
        - ((beta_1 - beta_2) * (beta_1 - 2*beta + beta_2) * eta_2) 
        + xi * omega * (beta + omega) * eta_c 
        - (beta_2 - beta)**2 * eta_3 
        - xi * omega * (beta + omega) 
        - n * (eta_s * rho_t + rho_c) * beta**2 
        + (kappa * n * eta_s * rho_t - kappa * n * rho_t - 2) * beta 
        + np.sqrt(
            (
                ((eta_s * beta - kappa * (eta_s - 1)) * rho_t + beta * rho_c)**2 * n**2 
                + (
                    (
                        2 * eta_s * (eta_c * xi * alpha - eta_3 * (alpha - 1)) * beta**2 
                        + (
                            ((4 * omega * (eta_c - 1) * alpha + (-2 * kappa - 2 * omega) * eta_c + 2 * omega) * xi 
                             - 4 * (alpha - 1) * common_term) * eta_s 
                            + 2 * eta_c * xi * kappa
                        ) * beta 
                        + (
                            2 * omega * (alpha * omega - kappa - omega) * (eta_c - 1) * xi 
                            + 2 * (alpha - 1) * beta_squared_terms
                        ) * eta_s 
                        + 2 * xi * kappa * omega * (eta_c - 1)
                    ) * rho_t 
                    - 2 * rho_c * (
                        (eta_c * (alpha - 1) * xi - alpha * eta_3) * beta**2 
                        + (
                            2 * (alpha - 0.5) * omega * (eta_c - 1) * xi 
                            - 2 * common_term * alpha
                        ) * beta 
                        + alpha * (
                            omega**2 * (eta_c - 1) * xi 
                            + beta_squared_terms
                        )
                    )
                ) * n 
                - xi * (
                    -eta_c * eta_3 * beta**2 
                    - 2 * eta_c * common_term * beta 
                    + omega**2 * (eta_c - 1) * xi 
                    + beta_squared_terms * eta_c
                ) * beta**2
            )
        ) + 1
    )
    
    denominator_k = (
        ((-1 + beta_1) * (beta_1 - 2*beta + 1) * eta_1) 
        - ((beta_1 - beta_2) * (beta_1 - 2*beta + beta_2) * eta_2) 
        + xi * (beta + omega)**2 * eta_c 
        - (beta_2 - beta)**2 * eta_3 
        + (-2 * beta * omega - omega**2) * xi 
        - 2 * beta 
        + 1
    )
    
    k = numerator_k / denominator_k
    
    # Compute M
    M_numerator = (
        (
            ((-2 * xi * eta_c + 2 * eta_3) * beta**3 
            + ((3 * eta_1 - 3 * eta_2) * beta_1 + (3 * eta_2 - 3 * eta_3) * beta_2 - 3 * eta_1 + 3 + (-3 * eta_c + 3) * omega * xi) * beta**2 
            + (-eta_1 + eta_2) * beta_1**3 + (-eta_2 + eta_3) * beta_2**3 + eta_1 - 1 + (eta_c - 1) * omega**3 * xi) * k**3 
            + (((-6 * eta_s * rho_t - 6 * rho_c) * n - 6 * eta_3) * beta**3 
            + (6 * rho_t * kappa * (eta_s - 1) * n + (-9 * eta_1 + 9 * eta_2) * beta_1 + (-9 * eta_2 + 9 * eta_3) * beta_2 + 9 * eta_1 - 9 + (3 * eta_c - 3) * omega * xi) * beta**2 
            + (3 * eta_1 - 3 * eta_2) * beta_1**3 + (3 * eta_2 - 3 * eta_3) * beta_2**3 - 3 * eta_1 + 3 + (-3 * eta_c + 3) * omega**3 * xi) * k**2 
            + (((12 * eta_s * rho_t * alpha - 12 * rho_c * (alpha - 1)) * n + 6 * eta_3) * beta**3 
            + (-6 * rho_t * kappa * (alpha + 1) * (eta_s - 1) * n + (9 * eta_1 - 9 * eta_2) * beta_1 + (9 * eta_2 - 9 * eta_3) * beta_2 - 9 * eta_1 + 9) * beta**2 
            + (-3 * eta_1 + 3 * eta_2) * beta_1**3 + (-3 * eta_2 + 3 * eta_3) * beta_2**3 + 3 * eta_1 - 3 + (3 * eta_c - 3) * omega**3 * xi) * k 
            + ((-6 * eta_s * rho_t * alpha**2 - 6 * rho_c * (alpha - 1)**2) * n - 2 * eta_3) * beta**3 
            + (6 * rho_t * alpha * kappa * (eta_s - 1) * n + (-3 * eta_1 + 3 * eta_2) * beta_1 + (-3 * eta_2 + 3 * eta_3) * beta_2 + 3 * eta_1 - 3) * beta**2 
            + (eta_1 - eta_2) * beta_1**3 + (eta_2 - eta_3) * beta_2**3 - eta_1 + 1 + (-eta_c + 1) * omega**3 * xi
        ) 
        * h**2 
        * epsilon_cr
    )
    
    M_denominator = 6 * (k - 1) * beta**2
    M = E * b * M_numerator / M_denominator
    
    return k, M


def zone422(beta_z4, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):

        beta = beta_z4
        
        k = ((-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1 - (beta_1 - beta_2) * (beta_1 - 2 * beta + beta_2) * eta_2 + 
                   xi * omega * (beta + omega) * eta_c - (beta_2 - beta)**2 * eta_3 - xi * omega * (beta + omega) - 
                   n * (eta_s * rho_t + rho_c) * beta**2 + (kappa * n * eta_s * rho_t - kappa * n * rho_t - 2) * beta + 
                   np.sqrt((((eta_s * beta - kappa * (eta_s - 1)) * rho_t + beta * rho_c)**2 * n**2 + 
                   ((2 * eta_s * (eta_c * xi * alpha - eta_3 * (alpha - 1)) * beta**2 + 
                   (((4 * omega * (eta_c - 1) * alpha + (-2 * kappa - 2 * omega) * eta_c + 2 * omega) * xi - 
                   4 * (alpha - 1) * ((eta_1 - eta_2) * beta_1 + (eta_2 - eta_3) * beta_2 - eta_1 + 1)) * eta_s + 
                   2 * eta_c * xi * kappa) * beta + 
                   (2 * omega * (alpha * omega - kappa - omega) * (eta_c - 1) * xi + 
                   2 * (alpha - 1) * ((eta_1 - eta_2) * beta_1**2 + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1)) * eta_s + 
                   2 * xi * kappa * omega * (eta_c - 1)) * rho_t - 
                   2 * rho_c * ((eta_c * (alpha - 1) * xi - alpha * eta_3) * beta**2 + 
                   (2 * (alpha - 0.5) * omega * (eta_c - 1) * xi - 
                   2 * ((eta_1 - eta_2) * beta_1 + (eta_2 - eta_3) * beta_2 - eta_1 + 1) * alpha) * beta + 
                   alpha * (omega**2 * (eta_c - 1) * xi + (eta_1 - eta_2) * beta_1**2 + 
                   (eta_2 - eta_3) * beta_2**2 - eta_1 + 1))) * n - 
                   xi * (-eta_c * eta_3 * beta**2 - 2 * eta_c * ((eta_1 - eta_2) * beta_1 + (eta_2 - eta_3) * beta_2 - 
                   eta_1 + 1) * beta + omega**2 * (eta_c - 1) * xi + 
                   ((eta_1 - eta_2) * beta_1**2 + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1) * eta_c)) * beta**2) + 1) / ((-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1 - (beta_1 - beta_2) * (beta_1 - 2 * beta + beta_2) * eta_2 + 
                    xi * (beta + omega)**2 * eta_c - (beta_2 - beta)**2 * eta_3 + (-2 * beta * omega - omega**2) * xi - 2 * beta + 1)

        M = (E * b * (((-2 * xi * eta_c + 2 * eta_3) * beta**3 + 
          ((3 * eta_1 - 3 * eta_2) * beta_1 + (3 * eta_2 - 3 * eta_3) * beta_2 - 3 * eta_1 + 3 + 
          (-3 * eta_c + 3) * omega * xi) * beta**2 + 
          (-eta_1 + eta_2) * beta_1**3 + (-eta_2 + eta_3) * beta_2**3 + eta_1 - 1 + (eta_c - 1) * omega**3 * xi) * k**3 + 
          (((-6 * eta_s * rho_t - 6 * rho_c) * n - 6 * eta_3) * beta**3 + 
          (6 * rho_t * kappa * (eta_s - 1) * n + (-9 * eta_1 + 9 * eta_2) * beta_1 + (9 * eta_2 - 9 * eta_3) * beta_2 + 
          9 * eta_1 - 9 + (3 * eta_c - 3) * omega * xi) * beta**2 + 
          (3 * eta_1 - 3 * eta_2) * beta_1**3 + (3 * eta_2 - 3 * eta_3) * beta_2**3 - 3 * eta_1 + 3 + 
          (-3 * eta_c + 3) * omega**3 * xi) * k**2 + 
          (((12 * eta_s * rho_t * alpha - 12 * rho_c * (alpha - 1)) * n + 6 * eta_3) * beta**3 + 
          (-6 * rho_t * kappa * (alpha + 1) * (eta_s - 1) * n + (9 * eta_1 - 9 * eta_2) * beta_1 + 
          (9 * eta_2 - 9 * eta_3) * beta_2 - 9 * eta_1 + 9) * beta**2 + 
          (-3 * eta_1 + 3 * eta_2) * beta_1**3 + (-3 * eta_2 + 3 * eta_3) * beta_2**3 + 3 * eta_1 - 3 + 
          (3 * eta_c - 3) * omega**3 * xi) * k + 
          ((-6 * eta_s * rho_t * alpha**2 - 6 * rho_c * (alpha - 1)**2) * n - 2 * eta_3) * beta**3 + 
          (6 * rho_t * alpha * kappa * (eta_s - 1) * n + (-3 * eta_1 + 3 * eta_2) * beta_1 + 
          (-3 * eta_2 + 3 * eta_3) * beta_2 + 3 * eta_1 - 3) * beta**2 + 
          (eta_1 - eta_2) * beta_1**3 + (eta_2 - eta_3) * beta_2**3 - eta_1 + 1 + 
          (-eta_c + 1) * omega**3 * xi) * h**2 * epsilon_cr / (6 * (k - 1) * beta**2))
        
        return k, M
'''

def zone4222(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    k = np.real(((-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1 - (beta_1 - beta_2) * (beta_1 - 2 * beta + beta_2) * eta_2 + xi * omega * (beta + omega) * eta_c - xi * omega**2 - beta * omega * xi - (beta_2 - beta)**2 * eta_3 - n * eta_s * (rho_c + rho_t) * beta**2 + (-2 - kappa * (rho_c - rho_t) * (eta_s - 1) * n) * beta + np.sqrt(beta**2 * ((((beta - kappa) * rho_t + rho_c * (beta + kappa)) * eta_s - (rho_c - rho_t) * kappa)**2 * n**2 + ((((2 * eta_c * xi * alpha - 2 * eta_3 * (alpha - 1)) * beta**2 + ((4 * omega * (eta_c - 1) * alpha + (-2 * kappa - 2 * omega) * eta_c + 2 * omega) * xi - 4 * ((eta_1 - eta_2) * beta_1 + (eta_2 - eta_3) * beta_2 - eta_1 + 1) * (alpha - 1)) * beta + 2 * omega * (alpha * omega - kappa - omega) * (eta_c - 1) * xi + 2 * ((eta_1 - eta_2) * beta_1**2 + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1) * (alpha - 1)) * rho_t - 2 * ((eta_c * (alpha - 1) * xi - alpha * eta_3) * beta**2 + ((2 * omega * (eta_c - 1) * alpha + (-kappa - omega) * eta_c + omega) * xi - 2 * ((eta_1 - eta_2) * beta_1 + (eta_2 - eta_3) * beta_2 - eta_1 + 1) * alpha) * beta + omega * (alpha * omega - kappa) * (eta_c - 1) * xi + ((eta_1 - eta_2) * beta_1**2 + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1) * alpha) * rho_c) * eta_s - 2 * (rho_c - rho_t) * (eta_c * beta + omega * (eta_c - 1)) * kappa * xi) * n - (-eta_c * eta_3 * beta**2 - 2 * ((eta_1 - eta_2) * beta_1 + (eta_2 - eta_3) * beta_2 - eta_1 + 1) * eta_c * beta + omega**2 * (eta_c - 1) * xi + ((eta_1 - eta_2) * beta_1**2 + (eta_2 - eta_3) * beta_2**2 - eta_1 + 1) * eta_c) * xi)) + 1) / ((-1 + beta_1) * (beta_1 - 2 * beta + 1) * eta_1 - (beta_1 - beta_2) * (beta_1 - 2 * beta + beta_2) * eta_2 + xi * (beta + omega)**2 * eta_c - xi * omega**2 - 2 * beta * omega * xi - (beta_2 - beta)**2 * eta_3 - 2 * beta + 1))
    M = (h**2 * E * epsilon_cr * b * (((-2 * xi * eta_c + 2 * eta_3) * beta**3 +
          ((3 * eta_1 - 3 * eta_2) * beta_1 + (3 * eta_2 - 3 * eta_3) * beta_2 - 3 * eta_1 + 3 +
          (-3 * eta_c + 3) * omega * xi) * beta**2 + (-eta_1 + eta_2) * beta_1**3 +
          (-eta_2 + eta_3) * beta_2**3 + eta_1 - 1 + (eta_c - 1) * omega**3 * xi) * k**3 +
          ((-6 * n * eta_s * (rho_c + rho_t) - 6 * eta_3) * beta**3 +
          (-6 * kappa * (rho_c - rho_t) * (eta_s - 1) * n + (-9 * eta_1 + 9 * eta_2) * beta_1 +
          (-9 * eta_2 + 9 * eta_3) * beta_2 + 9 * eta_1 - 9 + (3 * eta_c - 3) * omega * xi) * beta**2 +
          (3 * eta_1 - 3 * eta_2) * beta_1**3 + (3 * eta_2 - 3 * eta_3) * beta_2**3 - 3 * eta_1 + 3 +
          (-3 * eta_c + 3) * omega**3 * xi) * k**2 +
          ((-12 * eta_s * ((alpha - 1) * rho_c - rho_t * alpha) * n + 6 * eta_3) * beta**3 +
          (-6 * (eta_s - 1) * ((alpha - 2) * rho_c + rho_t * (alpha + 1)) * kappa * n +
          (9 * eta_1 - 9 * eta_2) * beta_1 + (9 * eta_2 - 9 * eta_3) * beta_2 - 9 * eta_1 + 9) * beta**2 +
          (-3 * eta_1 + 3 * eta_2) * beta_1**3 + (-3 * eta_2 + 3 * eta_3) * beta_2**3 + 3 * eta_1 - 3 +
          (3 * eta_c - 3) * omega**3 * xi) * k +
          (-6 * eta_s * ((alpha - 1)**2 * rho_c + rho_t * alpha**2) * n - 2 * eta_3) * beta**3 +
          (6 * ((alpha - 1) * rho_c + rho_t * alpha) * (eta_s - 1) * kappa * n +
          (-3 * eta_1 + 3 * eta_2) * beta_1 + (-3 * eta_2 + 3 * eta_3) * beta_2 + 3 * eta_1 - 3) * beta**2 +
          (eta_1 - eta_2) * beta_1**3 + (eta_2 - eta_3) * beta_2**3 - eta_1 + 1 + (-eta_c + 1) * omega**3 * xi) / (6 * (k - 1) * beta**2))

    return k, M
'''

import numpy as np
import math
def init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    global A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11
    global B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13
    global C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11
    global D1, D2, D3, D4, D5, D6, D7, D8, D9, D10
    global E1, E2, E3, E4, F1, F2, F3, F4

    # Shared global variables
    A1 = rho_t + rho_c
    A2 = (-2 * xi + 2) * alpha + 2 * xi
    A3 = 2 * rho_t * (1 + (xi - 1) * alpha)
    A4 = xi - 1
    A5 = -1 + (-rho_c - rho_t) * n
    A6 = (xi / 3 - 1 / 3)
    A7 = 1 + (rho_t + rho_c) * n
    A8 = -1 + ((2 * alpha - 2) * rho_c - 2 * rho_t * alpha) * n
    A9 = 1 / 3 + ((alpha - 1) ** 2 * rho_c + rho_t * alpha ** 2) * n
    A10 = 3 * (np.ones_like(beta) - 1)
    A11 = -3 * E * b * beta * epsilon_cr * h ** 2

    # Zone-specific global variables
    B1 = (-(-1 + beta) ** 2 * eta_1 + beta ** 2 * xi - 2 * beta + 1)
    B2 = (n * rho_t + xi / 2)
    B3 = (n * rho_c ** 2 + (2 * n * rho_t + 2 * xi) * rho_c + n * rho_t ** 2)
    B4 = (4 * n * rho_t + 2 * xi)
    B5 = (-eta_1 + (-rho_c - rho_t) * n)
    B6 = (2 * eta_1 - 2)
    B7 = (xi - eta_1) * beta ** 2 + B6 * beta - eta_1 + 1
    B8 = (-eta_1 / 3 + xi / 3)
    B9 = (eta_1 + (rho_c + rho_t) * n)
    B10 = (-eta_1 + ((2 * alpha - 2) * rho_c - 2 * rho_t * alpha) * n)
    B11 = ((alpha ** 2 - 2 * alpha + 1) * rho_c + rho_t * alpha ** 2) * n
    B12 = (np.ones_like(beta) - 1) ** 3 * (eta_1 - 1)
    B13 = beta ** 2 * (np.ones_like(beta) - 1)

    C1 = -(-1 + beta) ** 2 * eta_1
    C2 = eta_s * rho_t + rho_c
    C3 = kappa * n * eta_s * rho_t - kappa * n * rho_t - 2
    C4 = (eta_s * rho_t + rho_c) * beta - rho_t * kappa * (eta_s - 1)
    C5 = (-alpha + 1) * eta_1 + xi * alpha
    C6 = -alpha * eta_1 + xi * (alpha - 1)
    C7 = ((4 * alpha - 4) * eta_1 - 2 * xi * kappa - 4 * alpha + 4)
    C8 = ((alpha - 1) * eta_s * rho_t - rho_c * alpha)
    C9 = beta ** 2 * eta_1 + (-2 * eta_1 + 2) * beta + eta_1 - 1
    C10 = (-(-1 + beta) ** 2 * eta_1 + beta ** 2 * xi - 2 * beta + 1)
    C11 = (eta_1 - xi)

    D1 = (xi * eta_c - eta_3)
    D2 = (eta_c - 1) * omega * xi
    D3 = (eta_1 - eta_2)
    D4 = (eta_2 - eta_3)
    D5 = beta ** 2 * A1 ** 2 * n ** 2
    D6 = (D3 * beta_1 ** 2 + D4 * beta_2 ** 2 - eta_1 + 1)
    D7 = D1 * beta ** 2 + (2 * D2 + (-2 * eta_1 + 2 * eta_2) * beta_1 + (-2 * eta_2 + 2 * eta_3) * beta_2 + 2 * eta_1 - 2) * beta + omega ** 2 * (eta_c - 1) * xi + D6
    D8 = ((2 * eta_c * xi * alpha - 2 * eta_3 * (alpha - 1)) * beta ** 2)
    D9 = ((-2 * kappa - 2 * omega) * eta_c + 2 * omega) * xi
    D10 = beta ** 2 * (eta_c * beta + omega * (eta_c - 1)) * xi


    # Additional globals for zones 322, 411, and 412
    E1 = (xi * eta_c - eta_3)
    E2 = (eta_c - 1) * omega * xi
    E3 = (eta_1 - eta_2)
    E4 = (eta_2 - eta_3)

    # Additional global variables for Zones 421, 422, and 4222
    F1 = (beta + omega) * xi
    F2 = (beta_1 - 2 * beta + 1)
    F3 = (beta_1 - beta_2)
    F4 = (beta_2 - beta)


def zone111(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(A1 ** 2 * n ** 2 + (A2 * rho_c + A3) * n + xi)
    k = (sqrt_term + A5) / A4

    # Calculate M
    M = A11 * (A6 * k ** 3 + A7 * k ** 2 + A8 * k + A9) / A10

    return k, M


def zone211(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        (-2 * n * alpha * B1 * (rho_c - rho_t) + 2 * B2 * (-1 + beta) ** 2 * eta_1 +
         n * B3 * beta ** 2 + B4 * beta - 2 * n * rho_t - xi) * beta ** 2
    )
    k = np.real((sqrt_term + B5 * beta ** 2 + B6 * beta - eta_1 + 1) / B7)

    # Calculate M
    M = -epsilon_cr * E * b * h ** 2 * (
        ((B8 * k ** 3 + B9 * k ** 2 + B10 * k + eta_1 / 3 + B11) * beta ** 3 +
         B12 * beta ** 2 / 2 - B12 / 6) / B13
    )

    return k, M


def zone212(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        C4 ** 2 * n ** 2 + (
            (2 * C5 * eta_s * rho_t - 2 * rho_c * C6) * beta ** 2 +
            (C7 * eta_s * rho_t + 2 * xi * kappa - 4 * alpha * rho_c * (eta_1 - 1)) * beta -
            2 * C8 * (eta_1 - 1)
        ) * n + xi * C9
    )
    k = np.real((C1 - n * C2 * beta ** 2 + C3 * beta + sqrt_term + 1) / C10)

    # Calculate M
    M = C11 * h ** 2 * E * epsilon_cr * (
        ((C1 * k ** 3 + B9 * k ** 2 + B10 * k + A6) * beta ** 3 -
         3 * (A8 * k ** 2 + B5 * k + B6) * beta ** 2 / 2 + B12 / 2) / B13
    )

    return k, M


def zone221(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        beta ** 2 * (beta ** 2 * A1 ** 2 * n ** 2 +
                     ((-2 * C5 * rho_c * xi + 2 * C5 * eta_1) * beta ** 2 +
                      C7 * beta - 2 * A8) * n - xi * (B3 - 2 * B4 * beta + C8))
    )
    k = np.real((sqrt_term + B5 * beta ** 2 + B6 * beta + 1) / B7)

    # Calculate M
    M = A11 * (
        ((C8 * k ** 3 + C6 * k ** 2 + B5 * k + A9) * beta ** 3 -
         (C4 * k ** 2 + A5 * k) * beta ** 2 + B12) / B13
    )

    return k, M

def zone222(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        D5 + ((D8 + (D9 - 4 * A8) * beta + D10) * rho_t - 2 * ((D1 * alpha - D3) * beta ** 2 +
            (2 * (alpha - 0.5) * omega * D2 - 2 * D6 * alpha) * beta + omega * (alpha * omega - kappa) * D2 +
            D6 * alpha) * rho_c) * n - (D7 * xi))
    k = np.real((sqrt_term + B5 * beta ** 2 + B6 * beta + 1) / D7)

    # Calculate M
    M = A11 * (
        ((D1 * k ** 3 + D2 * k ** 2 + D3 * k + A9) * beta ** 3 -
         (C4 * k ** 2 + A5 * k) * beta ** 2 + B12) / B13
    )

    return k, M


def zone311(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):

    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        D5 + (D8 + (D9 - 4 * A8) * beta + D10) * rho_t - 2 * ((D1 * alpha - D3) * beta ** 2 +
            (2 * (alpha - 0.5) * omega * D2 - 2 * D6 * alpha) * beta + omega * (alpha * omega - kappa) * D2 +
            D6 * alpha) * rho_c * n - (D7 * xi))
    k = np.real((sqrt_term + B5 * beta ** 2 + B6 * beta + 1) / D7)

    # Calculate M
    M = A11 * (
        ((D1 * k ** 3 + D2 * k ** 2 + D3 * k + A9) * beta ** 3 -
         (C4 * k ** 2 + A5 * k) * beta ** 2 + B12) / B13
    )

    return k, M

def zone312(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        D5 + (D8 + (D9 - 4 * A8) * beta + D10) * rho_t - 2 * ((D1 * alpha - D3) * beta ** 2 +
            (2 * (alpha - 0.5) * omega * D2 - 2 * D6 * alpha) * beta + omega * (alpha * omega - kappa) * D2 +
            D6 * alpha) * rho_c * n - (D7 * xi))
    k = np.real((sqrt_term + B5 * beta ** 2 + B6 * beta + 1) / D7)

    # Calculate M
    M = A11 * (
        ((D1 * k ** 3 + D2 * k ** 2 + D3 * k + A9) * beta ** 3 -
         (C4 * k ** 2 + A5 * k) * beta ** 2 + B12) / B13
    )

    return k, M

def zone321(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        beta ** 2 * (
            D5 + (D8 + D9 * beta + D10) * rho_t - 
            2 * ((E1 * alpha - E3) * beta ** 2 +
                 (D9 - E4 * alpha) * beta + 
                 E4 * alpha) * rho_c * n - 
            D7 * xi
        )
    )
    k = np.real((sqrt_term + B5 * beta ** 2 + B6 * beta + 1) / D7)

    # Calculate M
    M = A11 * (
        ((E1 * k ** 3 + E2 * k ** 2 + E3 * k + A9) * beta ** 3 -
         (C4 * k ** 2 + A5 * k) * beta ** 2 + B12) / B13
    )

    return k, M

def zone322(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        D5 + (D8 + (D9 - 4 * A8) * beta + D10) * rho_t -
        2 * ((D1 * alpha - E3) * beta ** 2 +
             (2 * (alpha - 0.5) * omega * D2 - 2 * D6 * alpha) * beta +
             omega * (alpha * omega - kappa) * D2 + D6 * alpha) * rho_c * n -
        (D7 * xi)
    )
    k = np.real((sqrt_term + B5 * beta ** 2 + B6 * beta + 1) / D7)

    # Calculate M
    M = A11 * (
        ((D1 * k ** 3 + D2 * k ** 2 + D3 * k + A9) * beta ** 3 -
         (C4 * k ** 2 + A5 * k) * beta ** 2 + B12) / B13
    )

    return k, M


def zone411(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        D5 + (D8 + (D9 - 4 * A8) * beta + D10) * rho_t -
        2 * ((D1 * alpha - E3) * beta ** 2 +
             (2 * (alpha - 0.5) * omega * D2 - 2 * D6 * alpha) * beta +
             omega * (alpha * omega - kappa) * D2 + D6 * alpha) * rho_c * n -
        (D7 * xi)
    )
    k = np.real((sqrt_term + B5 * beta ** 2 + B6 * beta + 1) / D7)

    # Calculate M
    M = A11 * (
        ((D1 * k ** 3 + D2 * k ** 2 + D3 * k + A9) * beta ** 3 -
         (C4 * k ** 2 + A5 * k) * beta ** 2 + B12) / B13
    )

    return k, M


def zone412(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        D5 + (D8 + (D9 - 4 * A8) * beta + D10) * rho_t -
        2 * ((D1 * alpha - E3) * beta ** 2 +
             (2 * (alpha - 0.5) * omega * D2 - 2 * D6 * alpha) * beta +
             omega * (alpha * omega - kappa) * D2 + D6 * alpha) * rho_c * n -
        (D7 * xi)
    )
    k = np.real((sqrt_term + B5 * beta ** 2 + B6 * beta + 1) / D7)

    # Calculate M
    M = A11 * (
        ((D1 * k ** 3 + D2 * k ** 2 + D3 * k + A9) * beta ** 3 -
         (C4 * k ** 2 + A5 * k) * beta ** 2 + B12) / B13
    )

    return k, M


def zone421(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        beta ** 2 * (
            D5 + (D8 + (D9 - 4 * A8) * beta + D10) * rho_t -
            2 * ((D1 * alpha - D3) * beta ** 2 +
                 (2 * (alpha - 0.5) * omega * D2 - 2 * D6 * alpha) * beta +
                 omega * (alpha * omega - kappa) * D2 + D6 * alpha) * rho_c * n -
            (D7 * xi)
        )
    )
    k = np.real((sqrt_term + B5 * beta ** 2 + B6 * beta + 1) / D7)

    # Calculate M
    M = A11 * (
        ((D1 * k ** 3 + D2 * k ** 2 + D3 * k + A9) * beta ** 3 -
         (C4 * k ** 2 + A5 * k) * beta ** 2 + B12) / B13
    )

    return k, M

def zone422(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        beta ** 2 * (
            D5 + (D8 + (D9 - 4 * A8) * beta + D10) * rho_t -
            2 * ((D1 * alpha - E3) * beta ** 2 +
                 (2 * (alpha - 0.5) * omega * D2 - 2 * D6 * alpha) * beta +
                 omega * (alpha * omega - kappa) * D2 + D6 * alpha) * rho_c * n -
            (D7 * xi)
        )
    )
    k = np.real((sqrt_term + B5 * beta ** 2 + B6 * beta + 1) / D7)

    # Calculate M
    M = A11 * (
        ((D1 * k ** 3 + D2 * k ** 2 + D3 * k + A9) * beta ** 3 -
         (C4 * k ** 2 + A5 * k) * beta ** 2 + B12) / B13
    )

    return k, M

def zone4222(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t):
    init_globals(beta, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)

    # Calculate k
    sqrt_term = np.sqrt(
        D5 + (D8 + (D9 - 4 * A8) * beta + D10) * rho_t -
        2 * ((D1 * alpha - E3) * beta ** 2 +
             (2 * (alpha - 0.5) * omega * D2 - 2 * D6 * alpha) * beta +
             omega * (alpha * omega - kappa) * D2 + D6 * alpha) * rho_c * n -
        (D7 * xi)
    )
    k = np.real((sqrt_term + B5 * beta ** 2 + B6 * beta + 1) / D7)

    # Calculate M
    M = A11 * (
        ((D1 * k ** 3 + D2 * k ** 2 + D3 * k + A9) * beta ** 3 -
         (C4 * k ** 2 + A5 * k) * beta ** 2 + B12) / B13
    )

    return k, M

'''