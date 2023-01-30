# PCRF-Info.Science-2021-An Efficient and High-Quality Pansharpening Model Based on Conditional Random Fields

For reduced scale experiments (simulation), please run *_RS.m, for full scale (real) experiments, please run *_FS.m.

If you wish to cite this paper, please refer to:

    @article{yangEfficientHighqualityPansharpening2021,
      title = {An efficient and high-quality pansharpening model based on conditional random fields},
      volume = {553},
      issn = {00200255},
      url = {https://linkinghub.elsevier.com/retrieve/pii/S0020025520311403},
      doi = {10.1016/j.ins.2020.11.046},
      abstract = {Pansharpening fuses a low spatial resolution multi-spectral ({MS}) image with the corresponding panchromatic ({PAN}) image to obtain a high spatial resolution {MS} ({HRMS}) image. Traditional fusion methods may easily cause a spectral or spatial distortion when injecting details into an {MS} image. To preserve the spectral and spatial information, an efﬁcient pansharpening model based on conditional random ﬁelds ({CRFs}) is proposed. With this model, a state feature function is designed to force the {HRMS} image ﬁltered using a blur function to be consistent with the up-sampled {MS} image and retain the spectral ﬁdelity. To obtain a proper blur function, a new ﬁlter-acquisition method is proposed for the uniﬁed {CRF}-based model. Meanwhile, a transition feature function is deﬁned to enable the transition of {HRMS} pixels to follow the gradient of a {PAN} image and ensure the sharpness of the fused image. Considering the characteristics of the gradient domain, a total variation regularization is designed to make the gradient of the {HRMS} image sparse. Finally, the augmented Lagrangian function of the model is solved by employing an alternating direction method of the multipliers. Experiment results indicate that, compared with previous state-of-theart pansharpening methods, the proposed method can achieve the best fusion results with high computational efﬁciency.},
      pages = {1--18},
      journaltitle = {Inf. Sci.},
      author = {Yang, Yong and Lu, Hangyuan and Huang, Shuying and Fang, Yuming and Tu, Wei},
      urldate = {2021-01-03},
      date = {2021-04},
      langid = {english},	
}
