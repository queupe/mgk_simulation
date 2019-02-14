// Organize to calculates queue MMK
#include <algorithm>
#include <iostream>
#include <list>
#include <fstream>
#include <iomanip> // setprecision
#include <chrono>
#include <math.h>
#include <stdio.h>

#define _DIR_INPUT_FILES "./input_data/"
#define _DIR_OUTPUT_FILES "./output/"

using namespace std;

class Job {
        private:
                unsigned int id, job_type;
                long double time_input;

        public:
                Job(){id = 0; job_type = 0; time_input = 0.0;}
                Job(unsigned int my_id, long double my_time_input, unsigned int my_type){
                        id = my_id;
                        time_input = my_time_input;
                        job_type = my_type;
                }

                void setValues(unsigned int my_id, long double my_time_input, unsigned int my_type) {
                        id = my_id;
                        time_input = my_time_input;
                        job_type = my_type;
                }
                unsigned int getId()        { return id;         }
                unsigned int getType()      { return job_type;   }
                long double getTime_input() { return time_input; }


};

class JobFull: public Job {
        private:
                unsigned int processor;
                long double time_process, time_out;

        public:
                JobFull(unsigned int my_id, long double my_time_input, unsigned int my_type): Job(my_id, my_time_input, my_type){ }


                void setProcessor (unsigned int process_id)       { processor    = process_id;     }
                void setTimeProcess (long double my_time_process) { time_process = my_time_process;}
                void setTimeOut (long double my_time_out)         { time_out     = my_time_out;    }

                unsigned int getProcessor ()   { return processor;    }
                long double getTimeProcess ()  { return time_process; }
                long double getTimeOut ()      { return time_out;     }


};

class Processor {
    private:
        std::list<JobFull> jobs;
        unsigned int number;
        unsigned int occuped;
        long double *times;
        long double *weight;
        std::size_t weight_sz;
        std::size_t job_sz;
        long double moment1_tt, moment2_tt, total_tt;
        long double mean_tt, var_tt, std_tt;

        unsigned int insert_time(long double proc_time);
        unsigned int remove_old (long double inp_time);

    public:
        Processor()
        {
            number     = 0;
            occuped   = 0;
            weight     = NULL;
            times      = NULL;
            weight_sz  = 0;
            job_sz     = 0;
            moment1_tt = 0.0;
            moment2_tt = 0.0;
            total_tt   = 0.0;
            mean_tt    = 0.0;
            var_tt     = 0.0;
            std_tt     = 0.0;
        }

        Processor(unsigned int number_of_servers, long double *my_weight, std::size_t my_weight_sz)
        {
            number     = number_of_servers;
            occuped   = 0;
            weight     = my_weight;
            weight_sz  = my_weight_sz;
            times      = new long double[number_of_servers];
            job_sz     = 0;
            moment1_tt = 0.0;
            moment2_tt = 0.0;
            total_tt   = 0.0;
            mean_tt    = 0.0;
            var_tt     = 0.0;
            std_tt     = 0.0;
        }

        void setValues(unsigned int number_of_servers, long double *my_weight, std::size_t my_weight_sz)
        {
          number     = number_of_servers;
          occuped   = 0;
          weight     = my_weight;
          weight_sz  = my_weight_sz;
          times      = new long double[number_of_servers];
          job_sz     = 0;
          moment1_tt = 0.0;
          moment2_tt = 0.0;
          total_tt   = 0.0;
          mean_tt    = 0.0;
          var_tt     = 0.0;
          std_tt     = 0.0;
        }

        JobFull update (JobFull job);


        long double getMean()
        {
            if (number  == 0)  { return 0.0;                      }
            if (mean_tt == 0.0){ mean_tt = moment1_tt / total_tt; }
            return mean_tt;
        }

        long double getVar()
        {
            if (number  == 0)   { return 0.0;                                           }
            if (mean_tt == 0.0) { mean_tt = moment1_tt / total_tt;                      }
            if (var_tt  == 0.0) { var_tt = (moment2_tt/total_tt) - (mean_tt * mean_tt); }
            return  var_tt;
        }

        long double getStd()
        {
            if (number  == 0)   { return 0.0;                                           }
            if (mean_tt == 0.0) { mean_tt = moment1_tt / total_tt;                      }
            if (var_tt  == 0.0) { var_tt = (moment2_tt/total_tt) - (mean_tt * mean_tt); }
            if (std_tt  == 0.0) { std_tt = sqrt(var_tt);                                }
            return std_tt;
        }
};

unsigned int Processor::insert_time(long double proc_time)
{
  unsigned int j=0;
  long double aux;
  Processor::times[j] = proc_time;
  j++;
  while (j < Processor::occuped && Processor::times[j-1] > Processor::times[j]){
        aux                   = Processor::times[j-1];
        Processor::times[j-1] = Processor::times[j];
        Processor::times[j]   = aux;
        j++;
      }
  return j - 1;  //Entry position
}
unsigned int Processor::remove_old (long double inp_time)
{
  unsigned int pos_initial, pos, i = 0;
  while (i < Processor::occuped && Processor::times[i] < inp_time){
    i++;
  }
  pos_initial = 0;
  for (pos = i; pos < Processor::occuped; pos++){
    Processor::times[pos_initial++] = Processor::times[pos];
  }
  if (pos_initial > 0){
    return pos_initial;
  }else {
    return Processor::occuped;
  }

}

JobFull Processor::update(JobFull job)
{

    // JobFull job_full = JobFull (job.getId(), job.getTime_input(), job.getType());
    long double first_out_time = 0.0;
    unsigned int type_protect = 0;
    unsigned int i,j;

    long double aux;
    std::size_t sz = jobs.size();
    std::list<JobFull>::iterator it1, it2;
    std::list<long double>::iterator it;

    //Protection for null inicialization
    if (weight == NULL) {
      job.setTimeOut(0);
      job.setTimeProcess(0);
      job.setProcessor(0);
      return job;
    }

    // When there are empty cores
    if (occuped < number)
    {
        // Setting the enter time in one core of the processor
        job.setTimeProcess(job.getTime_input());
        occuped++;
    }else {
        if (job.getTime_input() > times[0]) {
          job.setTimeProcess(job.getTime_input());
        }else {
          job.setTimeProcess(times[0]);
        }
    }

    // Setting the out time of job, and protecting the type size
    type_protect = job.getType();
    if ( type_protect > weight_sz) { type_protect = weight_sz; }
    job.setTimeOut( (*(weight + (type_protect - 1)) * number) + job.getTimeProcess() );


    Processor::insert_time(job.getTimeOut());
    Processor::occuped = Processor::remove_old (job.getTime_input());

    // Setting the core number
    job.setProcessor(occuped);


    if (occuped > number) {
      cout << "Erro de número de servidores servers("<< number <<"), lista("<< occuped <<") " << endl;
    }

    aux = job.getTimeOut() - job.getTime_input();
    moment1_tt += aux;
    moment2_tt += aux*aux;
    total_tt++;

    return job;
}




int read_mmk_entry( long double Es, long double El, long double alpha,
                    long double rho, int round, int samples,
                    unsigned int max_cores)
{
    int i, k;
    long unsigned int job_id, job_type, total_events=0;
    long double job_input_time;
    long double aux, var, mean, std, moment1=0.0, moment2=0.0;
    char str[1024], str1[1024];

    long double weights[]   = { Es, El };
    int weights_sz = 2;

    // Create a list containing integers
    std::list<JobFull> JobList;
    std::list<JobFull>::iterator itj;
    std::list<Processor> ProcessorList, ProcessListAux;
    std::list<Processor>::iterator itp;


    std::string input_filename, output_filename;
    std::string line_read;
    std::size_t found;
    std::string::size_type sz1, sz2;

    Processor proc;
    std::ifstream myfile;
    std::ofstream myofile;

    auto start = chrono::steady_clock::now();

    //  Insert the code that will be timed
    sprintf(str, "B_factor_%6.4Lf00_alpha_%4.2Lf00_rho_%4.2Lf00_round_%d.csv",
                  Es/El, alpha, rho, round+1);
    output_filename = std::string(str);
    str[10] = str[25] = str[36] = str[49] = '_';
    input_filename  = std::string(str);

    cout << "Input file:\'"<< _DIR_INPUT_FILES + input_filename << "\'"<< endl;

    // Create the list of processor, from 0 until max_cores-1
    for (i = 0; i < max_cores; i++){
      Processor proc = Processor(i+1, weights, weights_sz);
      ProcessorList.push_back(proc);
    }


    myfile.open(_DIR_INPUT_FILES + input_filename);
    if (myfile.is_open()) {
        i = 0;
        while  ( getline (myfile,line_read) && (i < samples || samples == 0) )
        {
                // File structure to read
		//       50.000000,    50155.704089,        1.000000
                //012345678901234567890123456789012345678901234567890
                //          1         2         3         4
                //cout << line_read << '\n';

                job_id = std::stoul(line_read, &sz1);
                found = line_read.find(',');
                job_input_time = std::stod (line_read.substr(found+1), &sz2);
                found = line_read.find(',', found+sz2);
                //cout << "sz=" << found+1 << " substring=" << line_read.substr(found+1) << endl;
                job_type = std::stoul(line_read.substr(found+1), &sz1);

                //
                //k = 0;
                for (itp=ProcessorList.begin(); itp != ProcessorList.end(); itp++){
                        JobFull job_aux1 = JobFull(job_id, job_input_time, job_type);
                        JobList.push_back(itp->update(job_aux1));
                }
                i++;
                //if (i % 1000 == 0)
                  //cout << "read " << i << " lines"<< endl;
        }
        myfile.close();
    }



    // Testing the output
    i = 1;
    k = 1;
    for (itj= JobList.begin(); itj != JobList.end(); itj++) {

      if (itj->getId() != i)
      {
        //cout << endl;
        i++;
        k = 1;
      }

      std::ofstream fd;
      sprintf(str1,"./output/results_cpp_factor_%7.5Lf_servers_%d.csv", Es/El , k);
      fd.open(std::string(str1),ios::app);
      if (fd.is_open()) {
      //if (k ==3) {
      // Print results
        // << "cores:"<< k++
        sprintf(str1, "%6u, %20.8Lf, %25.8Lf, %20.8Lf, %2u, %2u, %20.8Lf, %20.8Lf, %16.2Lf",
              itj->getId()                                   ,
              itj->getTime_input()                           ,
              itj->getTimeProcess()                          ,
              itj->getTimeOut()                              ,
              itj->getType()                                 ,
              itj->getProcessor()                            ,
              itj->getTimeOut()     - itj->getTime_input()   ,
              itj->getTimeProcess() - itj->getTime_input()   ,
              itj->getTimeOut()     - itj->getTimeProcess()  );
        fd << std::string(str1) << endl;

        //cout << "\tId:"<< itj->getId()                        << ",\t"
        //     << itj->getTime_input()                          << ",\t"
        //     << itj->getTimeProcess()                         << ",\t"
        //     << itj->getTimeOut()                             << ",\t"
        //     << itj->getTimeOut()     - itj->getTime_input()  << ",\t"
        //     << itj->getTimeProcess() - itj->getTime_input()  << ",\t"
        //     << itj->getTimeOut()     - itj->getTimeProcess() << endl;
      }
      k++;

      i = itj->getId();
    }


    i = 1;
    myofile.open(_DIR_OUTPUT_FILES + output_filename);
    if (myofile.is_open())
    {
      for(itp=ProcessorList.begin(); itp != ProcessorList.end(); itp++){
          //cout << i << "\tMean:" << itp->getMean() << "\tVariance:"<< itp->getVar() << "\tStandard deviation:"<< itp->getStd() << endl;
          sprintf(str1, "%10d, %20.8Lf, %25.8Lf, %20.8Lf ",i, itp->getMean(), itp->getVar(), itp->getStd());
          myofile << std::string(str1) << endl;
          i++;
      }
    }
    auto end = chrono::steady_clock::now();
    // Store the time difference between start and end
    auto diff = end - start;

    aux = (long double) chrono::duration <double, milli> (diff).count();
    //var = aux - (aux % 1000);
    //std = aux % 1000;
    //cout << "Time elapsed: "<< var << " s" << std << " ms; " << aux << " ms" << endl;
    cout << "Time elapsed: "<< aux << " ms" << endl;
    //cout << chrono::duration <double, nano> (diff).count() << " ns" << endl;

   JobList.clear();
   ProcessorList.clear();
   return 0;
}

int main(int argc, char *argv[]){
    std::string filename;
    char str[1024];
    //ofstream myfile;
    //myfile.open (filename);
    //myfile << "Writing this to a file.\n";
    //myfile.close();

    long double Es                 = 54.13;
    long double El, El_vec[]       = {Es/0.0005, Es/0.005, Es/0.05};
    long double alpha, alpha_vec[] = {0.99     , 0.8     , 0.6    };
    long double rho, rho_vec[]     = {0.95     , 0.8     , 0.5    };

    int El_vec_sz           = 1;
    int alpha_vec_sz        = 1;
    int rho_vec_sz          = 1;
    int round, rounds_t     = 1;

    int samples             = 1000;

    int  iEl, ialpha, irho, l = 0, max_cores = 5;

    for (iEl = 0; iEl < El_vec_sz; iEl++)
    {
      for(ialpha = 0; ialpha < alpha_vec_sz; ialpha++)
      {
        for(irho = 0; irho < rho_vec_sz; irho++)
        {
          for(round = 0; round < rounds_t; round++)
          {
            l++;
            El    = El_vec[iEl];
            alpha = alpha_vec[ialpha];
            rho   = rho_vec[irho];
            //            B_factor_0_000500_alpha_0_9900_rho_0_9500_round_1.csv

            //printf("%d - %s", l, str);
            //cout << l << " - " << filename << endl;
            read_mmk_entry(Es, El, alpha, rho, round, samples, max_cores);
          }
        }
      }
    }




    //read_mmk_entry(filename);

    return 0;
}
