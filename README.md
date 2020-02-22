# Introduction
This repo. is a clone of the PRIDE-PPPAR software developed by PRIDE Lab of Wuhan University for better organizing the codes and data.

The aim of our project (for the course [GNSS Lab at ETH Zurich](http://www.vvz.ethz.ch/lerneinheitPre.do?semkez=2020S&lerneinheitId=135561&lang=en)) is to compare various freely available, non-commercial PPP software products (CNES PPP Wizard, PRIDE-PPPAR, GFZ GAMP etc.). Beisdes, we'd like to compare the performance of real-time PPP and RTK (by RTKLib) at different baseline lengths. 

## Intro. of  PRIDE-PPPAR

![pridelab.icon](https://github.com/PrideLab/PRIDE-PPPAR/blob/master/pride.png)

[PRIDE-PPPAR](http://pride.whu.edu.cn/newsDetails.shtml?newskindid=20190228093001384DTk8BHLcatWNl) is an open source software package aimed at post-processing of GPS data. 
The opensource software originates in Dr. Maorong Ge's efforts on PPP-AR and later developed and improved by Dr. Jianghui Geng. It is an open-source software package which is based on many GNSS professionals' collective work in GNSS Research Center, Wuhan University. 

[PRIDE Lab](http://pride.whu.edu.cn/indexone.shtml) at Wuhan University release this package with the hope to advance high-precision applications in geodetic and geophysical fields, such as crustal motion and troposphere sounding studies.

PRIDE-PPPAR are available for:

* Different mode:  
    * Static, PPP float solution
    * Static, PPP ambiguity resolution
    * Kinematic, PPP float solution
    * Kinematic, PPP ambiguity resolution
* High-rate GPS data
    * 1Hz, 5hz, 10Hz    

Notes: The phase clock/bias products, which are computed and released by Wuhan University in bias-SINEX format [download link](ftp://igs.gnsswhu.cn/pub/whu/phasebias), are required by PRIDE-PPPAR.

## License
***Copyright (C) 2019 by Wuhan University, All rights reserved.***

## Literature

If you find this work useful in your research, please consider cite:

```
@article{geng2019pride,
  title={PRIDE PPP-AR: an open-source software for GPS PPP ambiguity resolution},
  author={Geng, Jianghui and Chen, Xingyu and Pan, Yuanxin and Mao, Shuyin and Li, Chenghong and Zhou, Jinning and Zhang, Kunlun},
  journal={GPS Solutions},
  volume={23},
  number={4},
  pages={91},
  year={2019},
  publisher={Springer}
}
```

```
@inproceedings{geng2018phase,
  title={Phase bias product and open-source software for undifferenced ambiguity resolution at Wuhan University},
  author={Geng, Jianghui and Chen, X},
  booktitle={IGS workshop},
  year={2018}
}
```

```
@article{geng2019modified,
  title={A modified phase clock/bias model to improve PPP ambiguity resolution at Wuhan University},
  author={Geng, Jianghui and Chen, Xingyu and Pan, Yuanxin and Zhao, Qile},
  journal={Journal of Geodesy},
  volume={93},
  number={10},
  pages={2053--2067},
  year={2019},
  publisher={Springer}
}
```

## Related blogs
[(知乎)利用开源pride_pppar 进行高频单点定位计算教程](https://zhuanlan.zhihu.com/p/66662278)

[(知乎)PRIDE-PPPAR源码阅读小记](https://zhuanlan.zhihu.com/p/101144206)

### Inutuition
Let's make the realm of Geodesy more open.

---
