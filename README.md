# Storage-D

[![License](https://img.shields.io/badge/release-v1.0-green.svg)](https://github.com/DNAstorage-iSynBio/Storage-D)

**Storage-D** is a DNA data codec tool that integrates multiple codec algorithms.

## Background

We developed a new DNA data codec algorithm - Wukong [^1]. In order to facilitate the use of new algorithm researchers, we developed a general framework of codec, which can directly access the core algorithm into the whole process.

We also implement two algorithms proposed by Church *et al.*[^2] and Goldman *et al.*[^3] and introduce them into our framework. Likewise, users can integrate YYC [^4]  and DNA Fountain [^5] algorithm into the system, and their code is available on github. *([DNA-storage-YYC](https://github.com/ntpz870817/DNA-storage-YYC), [dna-fountain](https://github.com/TeamErlich/dna-fountain))*

Our system helps users to easily encode files into DNA sequences or decode DNA sequences into raw files using different algorithms. At the same time, users can also add their algorithms and use the common parts of the tool, such as splitting data, adding RS codes, etc.

## Environment

This tool is developed in python and uses some third-party modules. 

Details are as follows:

- python: 3.10.4
- tqdm: 4.64.0
- reedsolo: 1.5.4
- primer3-py: 0.6.1
- pandas: 1.4.3

At the same time, if you want to use the primer generation program, you need to install **[blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)** in the system.


## Directory tree

```html
├── StorageD
│   ├── __init__.py
│   ├── abstract_codec.py
│   ├── church.py
│   ├── churchDecode.py
│   ├── codec.py
│   ├── ecc.py
│   ├── getPrimerPair.py
│   ├── goldman.py
│   ├── goldmanDecode.py
│   ├── primer3Setting.py
│   ├── rules.py
│   ├── tools.py
│   └── wukong.py
├── README.md
└── tests
    ├── testFile
    │   └── dna.jpg
    ├── testResult
    │   ├── dna_church.fasta
    │   ├── dna_church_decode.jpg
    │   ├── dna_goldman.fasta
    │   ├── dna_goldman_decode.jpg
    │   ├── dna_wukong.fasta
    │   └── dna_wukong_decode.jpg
    └── test_codec.py
```

## Usage

You can refer to the usage in **tests**. 

Here is an example of *Wukong algorithm*.

#### Encode
```py
from StorageD.codec import WukongEncode
#encode
encode_worker = WukongEncode(input_file_path="testFile/dna.jpg", output_dir="testResult/", 
                            sequence_length=200,max_homopolymer=4,min_gc=0.4, max_gc=0.6,rule_num=1, 
                            rs_num=0, add_redundancy=True,add_primer=True, primer_length=20)
res_file = encode_worker.common_encode()
```
#### Decode
```py
from StorageD.codec import WukongDecode
#decode
decode_worker = WukongDecode(input_file_path='testResult/dna_wukong.fasta',
                            output_dir="testResult/", rule_num=1)
res_file = decode_worker.common_decode()
```

## Citing

[^1]: Huang X, Cui J, Qiang W, Ye J, Wang Y, Xie X, Li Y, Dai J. Storage-D: a user-friendly tool that enables practical and personalized data storage into DNA.

[^2]: [Church, G. M., Gao, Y. & Kosuri, S. Next-generation digital information storage in DNA. Science. 2012; 33 (6102):1628.](https://www.science.org/doi/10.1126/science.1226355)

[^3]: [Goldman, N. et al. Towards practical, high-capacity, low-maintenance information storage in synthesized DNA. Nature. 2013; 494(7435):77-80.](https://www.nature.com/articles/nature11875)

[^4]: [Erlich, Y. & Zielinski, D. DNA Fountain enables a robust and efficient storage architecture. Science. 2017; 355(6328):950-954.](https://www.science.org/doi/10.1126/science.aaj2038)

[^5]: [Ping Z., Chen S, Zhou G, Huang X, Zhu S, Zhang H, Lee H, Lan Z, Cui J, Chen T, Zhang W, Yang H, Xu X, Church G, Shen Y. Towards Practical and Robust DNA-Based Data Archiving Using ‘ Yin-Yang Codec ’ System. Nature Computational Science. 2022; 2, 234–242.](https://www.nature.com/articles/s43588-022-00231-2)