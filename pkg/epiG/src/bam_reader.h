/*
 * bamReader.hpp
 *
 *  Created on: Jan 25, 2014
 *      Author: martin
 */

#ifndef BAMREADER_HPP_
#define BAMREADER_HPP_

#include <vector>
#include <string>
#include "samtools/bam.h"
#include "samtools/sam.h"

//Aligned read

class aligned_read {

public:

	t_position position;
	t_length length;
	t_seq_bases bases;
	t_epsilon_quality epsilon;

	aligned_read(t_seq_bases const& bases, t_epsilon_quality const& epsilon, t_position position) :
		position(position), length(bases.n_elem), bases(bases), epsilon(epsilon) {}

	aligned_read(std::string const& bases, std::string const& quality, t_position position) :
		position(position), length(bases.size()), bases(create_bases_vector(bases)), epsilon(create_epsilon_vector(quality)) {}

    aligned_read() : position(0), length(0), bases(), epsilon() {}
};

class read_info {
public:

    t_position position;
    t_length length;

    read_info(t_position position, t_length length) : position(position), length(length) {}

};

// Call back fetch function
static int fetch_func(const bam1_t *b, void *data)
{
	std::vector<aligned_read> * reads = static_cast<std::vector<aligned_read>*>(data);

	int length = b->core.l_qseq;
	t_position pos = b->core.pos;

	uint8_t * s = bam1_seq(b);
	t_seq_bases bases(length);

	uint8_t * q = bam1_qual(b);
	t_epsilon_quality epsilon(length);

	for(int i = 0; i < length; ++i) {
        bases(i) = seq_base_to_int(static_cast<char>(bam1_seqi(s,i)));
		epsilon(i) = quality_to_epsilon(static_cast<char>(q[i]));
	}

	reads->push_back(aligned_read(bases, epsilon, pos));

    return 0;
}

// Call back count reads function
static int fetch_info_func(const bam1_t *b, void *data)
{
    std::vector<read_info> * infos = static_cast<std::vector<read_info>*>(data);

    int length = b->core.l_qseq;
    t_position pos = b->core.pos;


    infos->push_back(read_info(pos,length));

    return 0;
}


class bamReader {

	std::string file;
	std::string ref;
	t_position start;
	t_position end;
	std::vector<aligned_read> reads;
    std::vector<read_info> infos;

public:

	bamReader(std::string const& bam_file, std::string const& ref_name, t_position start_pos, t_position end_pos) :
		file(bam_file), ref(ref_name), start(start_pos), end(end_pos), reads() {}

	void fetch() {

		//Open bam file
		samfile_t *file_handle = samopen(file.c_str(), "rb", 0);

		if (file_handle == 0) {
			throw std::runtime_error("Fail to open BAM file.\n");
		}

		//Open index file
       bam_index_t *idx_handle = bam_index_load(file.c_str()); // load BAM index

       if (idx_handle == 0) {

    	   samclose(file_handle);

    	   throw std::runtime_error("BAM indexing file is not available.\n");
       }

       //Get ref id
       int ref_id;
       int tmp_s;
       int tmp_e;
       bam_parse_region(file_handle->header, ref.c_str(), &ref_id, &tmp_s, &tmp_e); // parse the region

       if (ref_id < 0) {

    	   bam_index_destroy(idx_handle);
    	   samclose(file_handle);

    	   throw std::runtime_error("Invalid region \n");
       }

       //Fetch reads
       bam_fetch(file_handle->x.bam, idx_handle, ref_id, start, end, &reads, fetch_func);

       bam_index_destroy(idx_handle);

       samclose(file_handle);
	}

    std::vector<aligned_read> const& get_reads() {
		return reads;
	}

    std::vector<read_info> const& get_infos() {
        return infos;
    }

    void fetch_reads_info() {

        //Open bam file
        samfile_t *file_handle = samopen(file.c_str(), "rb", 0);

        if (file_handle == 0) {
            throw std::runtime_error("Fail to open BAM file.\n");
        }

        //Open index file
       bam_index_t *idx_handle = bam_index_load(file.c_str()); // load BAM index

       if (idx_handle == 0) {

           samclose(file_handle);

           throw std::runtime_error("BAM indexing file is not available.\n");
       }

       //Get ref id
       int ref_id;
       int tmp_s;
       int tmp_e;
       bam_parse_region(file_handle->header, ref.c_str(), &ref_id, &tmp_s, &tmp_e); // parse the region

       if (ref_id < 0) {

           bam_index_destroy(idx_handle);
           samclose(file_handle);

           throw std::runtime_error("Invalid region \n");
       }

       //Count reads
       bam_fetch(file_handle->x.bam, idx_handle, ref_id, start, end, &infos, fetch_info_func);

       bam_index_destroy(idx_handle);

       samclose(file_handle);
    }
};


#endif /* BAMREADER_HPP_ */
