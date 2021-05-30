#include <iostream>
#include <fstream>
#include <zlib.h>

#include "kseq.h"
#include "Common.hpp"
#include "BUSData.h"

#include "bustools_countunmappedmolreads.h"
#include <map>
#include "BUSData.h"

KSEQ_INIT(gzFile, gzread);

void bustools_countunmappedmolreads(const Bustools_opt &opt) {
	BUSHeader h;
	size_t nr = 0;
	size_t N = 100000;
	uint32_t bclen = 0;
	uint32_t umilen = 0;
	uint32_t flag = 0;
	BUSData* p = new BUSData[N];

	std::vector<uint64_t> keys;
	keys.reserve(10^6);	
	std::vector<uint32_t> busMatches;
	busMatches.reserve(10^6);
	std::vector<uint32_t> fastqMatches;
	fastqMatches.reserve(10^6);
	
	uint64_t nMappedReads=0, nMatchedReads=0, nUnmatchedReads=0;

	//std::unordered_map<uint64_t, std::pair<int32_t,int32_t>> countsPerMol;//first is mapped counts from the bus file, second is unmapped counts
	//countsPerMol.reserve(10^6);

	//So, fill the map from the bus file first

	std::vector<BUSData> v;
	v.reserve(N);
	uint64_t current_bc = 0xFFFFFFFFFFFFFFFFULL;

	auto handleBarcode = [&](const std::vector<BUSData> &v) {
		if(v.empty()) {
		  return;
		}
		
		//create string barcode to be safe
		//std::string strBC;
		//binaryToString(v[0].barcode, bclen);

		double val = 0.0;
		size_t n = v.size();

		for (size_t i = 0; i < n; ) {
		  size_t j = i+1;
		  for (; j < n; j++) {
			if (v[i].UMI != v[j].UMI) {
			  break;
			}
		  }
		  
		  uint64_t combinedBCUMI = v[i].UMI + (v[i].barcode << (umilen * 2));
		  
		  std::string strUMI = binaryToString(v[i].UMI, umilen);
		  std::string strBC = binaryToString(v[i].barcode, bclen);
		  std::string strComb = binaryToString(combinedBCUMI, bclen + umilen);
		  //std::cout << umilen << '\t' << bclen << std::endl;
		  //std::cout << strBC << '\t' << strUMI << '\t'<< strComb << '\n'; 
		  
		  //auto combinedBCUMI = stringToBinary(strBC + strUMI, flag);
		  
		  
		  // v[i..j-1] share the same UMI
		  uint32_t counts = 0;
		  for (size_t k = i; k < j; k++) {
			counts += v[k].count;
		  }
		  
		  //countsPerMol.insert(std::make_pair(combinedBCUMI, std::pair<uint32_t,uint32_t>(counts, 0)));
		  keys.push_back(combinedBCUMI);
		  busMatches.push_back(counts);
		  fastqMatches.push_back(0);
		  nMappedReads += counts;
		  
		  i=j;
		}
	};

	for (const auto& infn : opt.files) { 
		std::streambuf *inbuf;
		std::ifstream inf;
		if (!opt.stream_in) {
		  inf.open(infn.c_str(), std::ios::binary);
		  inbuf = inf.rdbuf();
		} else {
		  inbuf = std::cin.rdbuf();
		}
		std::istream in(inbuf); 

		parseHeader(in, h);
		bclen = h.bclen;
		umilen = h.umilen;

		std::cout << "Reading BUS files\n";
		
		int rc = 0;
		while (true) {
		  in.read((char*)p, N*sizeof(BUSData));
		  size_t rc = in.gcount() / sizeof(BUSData);
		  nr += rc;
		  if (rc == 0) {
			break;
		  }
		  for (size_t i = 0; i < rc; i++) {
			if (p[i].barcode != current_bc) {                 
			  // output whatever is in v
			  if (!v.empty()) {
				  handleBarcode(v);
			  }
			  v.clear();
			  current_bc = p[i].barcode;
			}
			v.push_back(p[i]);

		  }            
		}
	
		if (!v.empty()) {
		  handleBarcode(v);
		}

		if (!opt.stream_in) {
		  inf.close();
		}
	}
	delete[] p; p = nullptr;

	uint64_t numFastqPairs = opt.fastq.size()/2;
	
	for (size_t iFastq = 0; iFastq < numFastqPairs; ++iFastq) {
		std::string fastq1Fn = opt.fastq[iFastq * 2];
		std::string fastq2Fn = opt.fastq[iFastq * 2 + 1];
		std::cout << "Reading fastq pair " << fastq1Fn << " " << fastq2Fn << "\n";
	  
		gzFile fp;  
		kseq_t *seq = nullptr;  
		int l1 = 0;  
//		int l2 = 0;
		auto fp1 = gzopen(fastq1Fn.c_str(), "r"); // STEP 2: open the file handler  
		auto seq1 = kseq_init(fp1); // STEP 3: initialize seq  
//		auto fp2 = gzopen(fastq2Fn.c_str(), "r"); // STEP 2: open the file handler  
//		auto seq2 = kseq_init(fp2); // STEP 3: initialize seq  

		while ((l1 = kseq_read(seq1)) >= 0) { // STEP 4: read sequence  
//			l2 = kseq_read(seq2);
	//		std::cout << "l1: " << l1 << " l2: " << l2 << "\n";
//			if (l2 < 1) {
//				std::cout << "Out of synch!\n";
//			}
			//std::cout << seq1->seq.s << std::endl;

			uint64_t bcumi = stringToBinary(seq1->seq.s, bclen + umilen, flag);
			//auto it = countsPerMol.find(bcumi);

			auto it = std::lower_bound(keys.begin(), keys.end(), bcumi);
			if ((it != keys.end()) && (*it == bcumi)) {
				fastqMatches[it - keys.begin()]++;
				nMatchedReads++;
			} else {
				nUnmatchedReads++;
			}
		}
		kseq_destroy(seq1); // STEP 5: destroy seq  
		gzclose(fp1); // STEP 6: close the file handler  		
//		kseq_destroy(seq2); // STEP 5: destroy seq  
//		gzclose(fp2); // STEP 6: close the file handler  		
	}

	std::cout << "Matched reads: " << nMatchedReads << " Unmatched reads: " << nUnmatchedReads << "\n";
	std::cout << "Writing output\n";

	
	//Now write the results to a text file
	std::ofstream of(opt.output);
	
	//for (auto p : countsPerMol) {
	for (size_t i = 0; i < keys.size(); ++i) {
		std::string totStr = binaryToString(keys[i], h.bclen + h.umilen);
		of << totStr.substr(0,bclen) << '\t' << totStr.substr(bclen) << '\t' << busMatches[i] << '\t' << fastqMatches[i] << '\n';
	}

	std::cout << "Matched reads: " << nMatchedReads << " Unmatched reads: " << nUnmatchedReads << "\n";
}

