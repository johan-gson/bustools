#include <iostream>
#include <fstream>
#include <zlib.h>

#include "kseq.h"
#include "Common.hpp"
#include "BUSData.h"

#include "bustools_extract.h"
#include <map>
#include "BUSData.h"

KSEQ_INIT(gzFile, gzread);

void bustools_extractforexon(const Bustools_opt &opt) {
	//read the molecules to save
	std:: ifstream ifs(opt.files[0]);
	//AAACCCAAGGCATCAG	GCTTATTAGGTT	4	Glutamatergic
	uint32_t flag = 0;
	
	struct FilterRec {
		std::string bc, umi;
		double counts;
		std::string cluster;
	};
	FilterRec rec;
	std::map<uint64_t, FilterRec> filterMap;
//	int i = 0;
	while (ifs >> rec.bc >> rec.umi >> rec.counts >> rec.cluster) {
		filterMap[stringToBinary(rec.bc + rec.umi, flag)] = rec;
//		if ( ++i < 20 ) {
//			std::cout << rec.bc << rec.umi << "\n";
//		}
	}
	
	std::ofstream of(opt.output);
  
	std::string fastq1Fn = opt.fastq[0];
	std::string fastq2Fn = opt.fastq[1];
  
    gzFile fp;  
    kseq_t *seq = nullptr;  
    int l1 = 0;  
	int l2 = 0;
    auto fp1 = gzopen(fastq1Fn.c_str(), "r"); // STEP 2: open the file handler  
    auto seq1 = kseq_init(fp1); // STEP 3: initialize seq  
    auto fp2 = gzopen(fastq2Fn.c_str(), "r"); // STEP 2: open the file handler  
    auto seq2 = kseq_init(fp2); // STEP 3: initialize seq  

    while ((l1 = kseq_read(seq1)) >= 0) { // STEP 4: read sequence  
		l2 = kseq_read(seq2);
//		std::cout << "l1: " << l1 << " l2: " << l2 << "\n";
		if (l2 < 1) {
			std::cout << "Out of synch!\n";
		}
		uint64_t bcumi = stringToBinary(seq1->seq.s, flag);
		auto it = filterMap.find(bcumi);
		if (it != filterMap.end()) {
//			std::cout << "Found: " << seq1->seq.s << "\n";
			of << it->second.bc << '\t' << it->second.umi << '\t' << it->second.counts << '\t' << it->second.cluster << '\t' << seq2->seq.s << '\n';
		} //else {
//			std::cout << "Not found: " << seq1->seq.s << "\n";
//		}
    }  
  
 /* 
  std::vector<gzFile> outFastq(opt.nFastqs);
  std::vector<gzFile> inFastq(opt.nFastqs);
  std::vector<kseq_t *> seq(opt.nFastqs, nullptr);
  uint32_t iRead = 0;
  size_t iFastq = 0;
  if (!open_fastqs(outFastq, inFastq, seq, opt, iFastq)) {
    std::cerr << "Error reading FASTQ " << opt.fastq[iFastq] << std::endl;
    goto end_extract;
  }

  while (true) {
    in.read((char *) p, N * sizeof(BUSData));
    size_t rc = in.gcount() / sizeof(BUSData);
    if (rc == 0) {
      break;
    }
    nr += rc;
    for (size_t i = 0; i < rc; ++i) {
      while (iRead < p[i].flags) {
        for (const auto &s : seq) {
          int err_kseq_read = kseq_read(s);
          if (err_kseq_read == -1) { // Reached EOF
            if (iFastq == opt.fastq.size()) { // Done with all files
              std::cerr << "Warning: number of reads in FASTQs was less than number of reads in BUS file" << std::endl;
              goto end_extract;
            } else {
              if (!open_fastqs(outFastq, inFastq, seq, opt, iFastq)) {
                std::cerr << "Error: cannot read FASTQ " << opt.fastq[iFastq] << std::endl;
                goto end_extract;
              }
            }
          } else if (err_kseq_read == -2) {
            std::cerr << "Error: truncated FASTQ" << std::endl;
            goto end_extract;
          }
        }
        ++iRead;
      }

      if (iRead > p[i].flags) {
        std::cerr << "BUS file not sorted by flag" << std::endl;
        goto end_extract;
      }

      for (int i = 0; i < opt.nFastqs; ++i) {
        int bufLen = 1; // Already have @ character in buffer
        
        memcpy(buf + bufLen, seq[i]->name.s, seq[i]->name.l);
        bufLen += seq[i]->name.l;
        
        memcpy(buf + bufLen, seq[i]->comment.s, seq[i]->comment.l);
        bufLen += seq[i]->comment.l;
        
        buf[bufLen++] = '\n';

        memcpy(buf + bufLen, seq[i]->seq.s, seq[i]->seq.l);
        bufLen += seq[i]->seq.l;
        
        buf[bufLen++] = '\n';
        buf[bufLen++] = '+';

        memcpy(buf + bufLen, seq[i]->name.s, seq[i]->name.l);
        bufLen += seq[i]->name.l;
        
        memcpy(buf + bufLen, seq[i]->comment.s, seq[i]->comment.l);
        bufLen += seq[i]->comment.l;
        
        buf[bufLen++] = '\n';

        memcpy(buf + bufLen, seq[i]->qual.s, seq[i]->qual.l);
        bufLen += seq[i]->qual.l;
        
        buf[bufLen++] = '\n';

        if (gzwrite(outFastq[i], buf, bufLen) != bufLen) {
          std::cerr << "Error writing to FASTQ" << std::endl;
          goto end_extract;
        }
      }
    }
  }

  std::cout << "Read in " << nr << " BUS records" << std::endl;

end_extract:
  delete[] p;
  delete[] buf;
  for (auto &elt : outFastq) {
    gzclose(elt);
  }
  for (auto &elt : inFastq) {
    gzclose(elt);
  }
  for (auto &elt : seq) {
    if (elt) {
      kseq_destroy(elt);
    }
  }
  
  */
  
  
}

