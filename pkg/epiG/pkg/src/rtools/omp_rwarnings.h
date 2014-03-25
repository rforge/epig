#ifndef OMP_RWARNINGS_H
#define OMP_RWARNINGS_H

class omp_rwarn {

private:

    std::vector<std::string> warn_msgs;

    void add_warning(std::string const& msg) {
 #ifdef EPIG_USE_OPENMP
 #pragma omp critical
 #endif
         {
             warn_msgs.push_back(msg);
         }
     }

public:

    omp_rwarn() {}

    ~omp_rwarn() {
        for(unsigned int i = 0; i < warn_msgs.size(); ++i) {
            Rf_warning(warn_msgs[i].c_str());
        }
    }

    void add(std::string const& msg) {
        add_warning(msg);
    }

};

#endif // OMP_RWARNINGS_H
