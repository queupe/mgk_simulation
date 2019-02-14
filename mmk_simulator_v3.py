import fileinput
import argparse
import math
import queue
import time
import sys
import re
import matplotlib.pyplot as plt
import numpy as np

from datetime import datetime
from pathlib import Path


def read_entry (filename = 'entry_jobs.csv'):
    entry_lst = list()
    with open(filename) as fd_in:
        for line in fd_in:
            entry_lst.append([int(float(line.strip(' ').split(',')[0])), float(line.strip(' ').split(',')[1]), int(float(line.strip(' ').rstrip('\n').split(',')[2]))])

    return entry_lst

def read_entry_max (filename = 'entry_jobs.csv', max = 100):
    entry_lst = list()
    count = 0
    with open(filename) as fd_in:
        for line in fd_in:
            count += 1
            entry_lst.append([int(float(line.strip(' ').split(',')[0])), float(line.strip(' ').split(',')[1]), int(float(line.strip(' ').rstrip('\n').split(',')[2]))])
            if count >= max:
                return entry_lst

    return entry_lst

def queue_mmk_typed (entry_lst, time_service_type, k_servers):
    job_hist = list() # id, time_enter, time_served, time_out, job_type, server
    serves_lst = list()

    for job in entry_lst:
        if len(serves_lst) < k_servers:
            job_id          = job[0]
            job_time_enter  = job[1]
            job_type        = job[2]

            job_time_server = job[1]
            job_time_out    = job_time_server + time_service_type[job[2]-1] * k_servers
            job_server      = len(serves_lst)

            job_hist.append([job_id, job_time_enter, job_time_server, job_time_out, job_type, job_server])
            serves_lst.append([job_id, job_time_enter, job_time_server, job_time_out, job_type, job_server])
        else:
            job_id          = job[0]
            job_time_enter  = job[1]
            job_type        = job[2]

            pos = np.argmin(np.array(serves_lst)[:,3])

            if serves_lst[pos][3] > job_time_enter:
                job_time_server = serves_lst[pos][3]
            else:
                job_time_server = job_time_enter

            job_time_out    = job_time_server + time_service_type[job[2]-1] * k_servers

            del serves_lst[pos]

            num_jobs = len(serves_lst)
            for i in range(num_jobs):
                pos = num_jobs - i -1
                if serves_lst[pos][3] < job_time_enter:
                    del serves_lst[pos]

            job_server = len(serves_lst)
            job_hist.append([job_id, job_time_enter, job_time_server, job_time_out, job_type, job_server])
            serves_lst.append([job_id, job_time_enter, job_time_server, job_time_out, job_type, job_server])

    #if (k_servers == 3):
    '''
    if True:
        with open("./output/results_python_factor_{:7.5f}_servers_{:2d}.csv".format(time_service_type[0]/time_service_type[1], k_servers),'a') as fd:
            for job_t in job_hist:
                fd.write('{:6d}, {:20.8f}, {:25.8f}, {:20.8f}, {:2d}, {:2d}, {:20.8f}, {:20.8f}, {:16.2f}\n'.format(job_t[0],
                    job_t[1]             ,
                    job_t[2]             ,
                    job_t[3]             ,
                    job_t[4]             ,
                    len(serves_lst)      ,
                    job_t[3] - job_t[1]  ,
                    job_t[2] - job_t[1]  ,
                    job_t[3] - job_t[2] ))
    '''
    return job_hist # id, time_enter, time_served, time_out, job_type, server

def queue_mmk_typed_with_error(entry_lst_vec, time_service_type, num_servers):

    job_servers_lst = list()
    rounds = 0
    for entry_lst in entry_lst_vec:
        h  = queue_mmk_typed(entry_lst, time_service_type, num_servers)
        h1 = np.array(h)
        rounds += 1

        # Values of: ....................Total time in sys ........ Queue time ..............Served time
        job_servers_lst.append([ np.mean(h1[:,3]-h1[:,1]), np.mean(h1[:,2]-h1[:,1]), np.mean(h1[:,3]-h1[:,2]),
                                 np.var (h1[:,3]-h1[:,1]), np.var (h1[:,2]-h1[:,1]), np.var (h1[:,3]-h1[:,2]),
                                 np.std (h1[:,3]-h1[:,1]), np.std (h1[:,2]-h1[:,1]), np.std (h1[:,3]-h1[:,2]),
                                 np.min (h1[:,3]-h1[:,1]), np.min (h1[:,2]-h1[:,1]), np.min (h1[:,3]-h1[:,2]),
                                 np.max (h1[:,3]-h1[:,1]), np.max (h1[:,2]-h1[:,1]), np.max (h1[:,3]-h1[:,2]) ])

    jb_array = np.mean(np.array(job_servers_lst), 0)
    error    = (1.96 / np.sqrt(rounds)) * np.std(np.array(job_servers_lst), 0)
    jb_array_error_p = jb_array + error
    jb_array_error_n = jb_array - error
    jb_array_min     = np.min(np.array(job_servers_lst), 0)
    jb_array_max     = np.max(np.array(job_servers_lst), 0)

    return (jb_array, jb_array_error_p, jb_array_error_n, jb_array_min, jb_array_max)



def simul_queue_servers (entry_lst, time_service_type, k_max_servers):

    job_servers_lst = list()
    for i in range (1, k_max_servers+1):
        h  = queue_mmk_typed(entry_lst, time_service_type, i)
        h1 = np.array(h)
        alpha = time_service_type[0]
        beta  = time_service_type[1]

        job_servers_lst.append([ np.mean(h1[:,3]-h1[:,1]), np.mean(h1[:,2]-h1[:,1]), np.mean(h1[:,3]-h1[:,2]),
                                 np.var (h1[:,3]-h1[:,1]), np.var (h1[:,2]-h1[:,1]), np.var (h1[:,3]-h1[:,2]),
                                 np.std (h1[:,3]-h1[:,1]), np.std (h1[:,2]-h1[:,1]), np.std (h1[:,3]-h1[:,2]),
                                 np.min (h1[:,3]-h1[:,1]), np.min (h1[:,2]-h1[:,1]), np.min (h1[:,3]-h1[:,2]),
                                 np.max (h1[:,3]-h1[:,1]), np.max (h1[:,2]-h1[:,1]), np.max (h1[:,3]-h1[:,2]) ])

    return job_servers_lst

def simul_queue_servers_with_erros (entry_lst_vec, time_service_type, k_max_servers):

    job_servers_lst_with_error = list()
    for num_servers in range (1, k_max_servers+1):
        jb  = queue_mmk_typed_with_error(entry_lst_vec, time_service_type, num_servers)
        job_servers_lst_with_error.append(jb)

    return job_servers_lst_with_error



def plotting_results(job_servers_lst_with_error, Es, El_factor, rho, dim):

    colors_vec = ['r', 'g', 'b']
    for d in dim:
        # Plotting the mean
        plt.fill_between(k[d], error_mu_p[d]   , error_mu_n[d]   , color=colors_vec[d%(len(colors_vec))], alpha=0.5)
        plt.fill_between(k[d], error_sigma_p[d], error_sigma_n[d], color=colors_vec[d%(len(colors_vec))], alpha=0.5)
        plt.semilogy(k[d], mu_t[d]   , label='E(T)'         , color=colors_vec[d%(len(colors_vec))])
        plt.semilogy(k[d], sigma_t[d], label='\\sigma(T)'   , color=colors_vec[d%(len(colors_vec))] , linestyle='--')

    plt.xlabel('number of servers (K)')
    plt.ylabel('response time')
    #plt.ylim((0.0,1.7))
    plt.legend(loc='upper left')
    fig_round_name ='img/fig_Bfactor_{:6.4f}_alpha_{:4.2f}_rho_{:4.2f}.png'.format(factor, alpha, rho)
    plt.savefig(fig_round_name, dpi=100)
    #plt.show()
    plt.clf()
    #copy(_WEB_LOCAL_ + fig_round_name, _WEB_REMOTE_+ fig_round_name ) # shutil


def main():

    alpha = 0.99
    factor= 0.0005
    rho = 0.50
    start_t = datetime.now()

    DIR_DATA = 'input_data/'

    Es = 54.13
    El_vec    = [Es/0.0005, Es/0.005, Es/0.05]
    alpha_vec = [0.6, 0.8, 0.99]
    rho_vec   = [0.95, 0.8, 0.5]
    rounds_t  = 5

    counter = 0

    num_servers = 100

    #only to tests -
    #Es = 54.13
    #El_vec    = [Es/0.0005]
    #alpha_vec = [0.99]
    #rho_vec   = [0.95]
    rounds_t  = 5

    num_samples = 100000
    #print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('serves','mean', 'err+', 'err-', 'sigma', 'err+', 'err-'))

    # Específic values:
    #Es = 54.13
    #----------->   El              alpha   rho
    #especific = (   (Es/0.0005,     0.8,    0.95),
    #                (Es/0.0005,     0.99,   0.95),
    #                (Es/0.0005,     0.99,   0.8 ),
    #                (Es/0.005,      0.99,   0.95))
    #num_samples = 500000

    for El in El_vec:
    #if True:
        for alpha in alpha_vec:
        #if True:
            for rho in rho_vec:
            #for espec in especific:
            #    (El, alpha, rho) = espec

                entry_lst_vec = list()
                for round in range(1,rounds_t+1):
                    input_file_entry = 'B_factor_{:8.6f}_alpha_{:6.4f}_rho_{:6.4f}_round_{:d}_csv'.format(Es/El, alpha, rho, round)
                    input_file_entry = re.sub('[.]','_',input_file_entry)
                    #print(input_file_entry)
                    entry_lst_vec.append(read_entry_max(DIR_DATA + input_file_entry, num_samples))

                #return 0

                time_server_type = [Es, El]
                results = simul_queue_servers_with_erros(entry_lst_vec, time_server_type, num_servers)
                k_aux = 0
                k_servers = range(1, len(results) + 1)

                file_round_name ='output/fig_Bfactor_{:6.4f}_alpha_{:4.2f}_rho_{:4.2f}.csv'.format(Es/El, alpha, rho)
                with open(file_round_name,'w') as fd_out:
                    for res in results:
                        k_aux +=1
                        #print('{:d}\t{}\t{}\t{}\t{}\t{}\t{}'.format(k_aux,res[0][0], res[1][0], res[2][0], res[0][6], res[1][6], res[2][6]))
                        fd_out.write('{:d},{},{},{},{},{},{}\n'.format(k_aux,res[0][0], res[1][0], res[2][0], res[0][6], res[1][6], res[2][6]))

                r = np.array(results)
                plt.fill_between(k_servers, r[:,1,0]   , r[:,2,0]            , color='b', alpha=0.3)
                plt.fill_between(k_servers, r[:,1,6]   , r[:,2,6]            , color='r', alpha=0.3)
                plt.semilogy    (k_servers, r[:,0,0]   , label='E(T)'        , color='b')
                plt.semilogy    (k_servers, r[:,0,6]   , label=r'$\sigma(T)$', color='r' , linestyle='--')
                plt.xlabel('number of servers (K)')
                plt.ylabel('response time')
                #plt.ylim((0.0,1.7))
                plt.title(r'$E_S={}; E_S/E_L={}; \alpha={}; \rho={}$'.format(Es, Es/El, alpha, rho))
                plt.legend(loc='upper left')
                fig_round_name ='img/fig_Bfactor_{:6.4f}_alpha_{:4.2f}_rho_{:4.2f}.png'.format(Es/El, alpha, rho)
                plt.savefig(fig_round_name, dpi=100)
                counter += 1
                print('{:2d} - Save file {}'.format(counter, fig_round_name))
                time_round = datetime.now() - start_t
                print("Finish this round with: {0}".format(time_round))
                start_t = datetime.now()
                #plt.show()
                plt.clf()

    return 0




if __name__ == '__main__':
	start_time = datetime.now()
	#print(sys.argv[1])
	results = 0
	results = main()
	time_elapsed = datetime.now() - start_time
	print("Time elapsed: {0}".format(time_elapsed))
	sys.exit(results)
