MMK Simulator

Usage:
        python3 mmk_simulator_v2.py

It expect read the input files in ./input_data with folow configuration:
        Es    = 54.13
        El    = [Es/0.0005, Es/0.005, Es/0.05]
        alpha = [0.6, 0.8, 0.99]
        rho   = [0.95, 0.8, 0.5]
        round = [1:5]

        'B_factor_{:8.6f}_alpha_{:6.4f}_rho_{:6.4f}_round_{:d}_csv'.format(Es/El, alpha, rho, round)

        replacing '.' with '_'

        The total files are:
        3 x 3 x 3 x 5 = 135 files

The format of input files is:
        #<Job ID>, #<input time>, #<Job type>
        
        example:
        1.000000,     1169.211830,        2.000000
        2.000000,     1892.379238,        1.000000
        3.000000,     2227.429930,        1.000000
        4.000000,     2602.631844,        2.000000
        5.000000,     3112.076599,        1.000000

The output is the:
        standard output:
                - Information about the made graphical files
                - and the time taken to generate the files 
        files in ".img/"
                - The graphical files with the same information of input files, 
                without replacements:
                
                'fig_Bfactor_{:6.4f}_alpha_{:4.2f}_rho_{:4.2f}.png'.format(Es/El, alpha, rho)
                
The file "gera_dados.m" is a matlab script to generate the input files.
        It save all files, in the format expected by "MMK Simulator", on the folder "./input_data/";
        This folder must exist.
     
