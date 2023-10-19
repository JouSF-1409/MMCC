# MMcc

## EN intro

modified from Jiang,2012. A phase weighting multichannel cross-correlation strategy to pick relative arrival time.

MMCC_test.jupyter contains full function, but without a proper management.

original formula is 

$$ 
C_{\mathrm{mmmcc}}(\tau)=\frac{\sum_{k=1}^{N} S_{j}(t+k \cdot \delta t-\tau) \cdot S_{k}(t+k \cdot \delta t) \cdot\left|\cos \left(\frac{\phi_{j}(t+k \cdot \delta t-\tau)-\phi_{k}(t+k \cdot \delta t)}{2}\right)\right|^{v}}{\sqrt{\sum_{k=1}^{N}\left|S_{j}(t+k \cdot \delta t-\tau)\right|^{2} \cdot \sum_{k=1}^{N}\left|S_{k}(t+k \cdot \delta t)\right|^{2}}}
$$

setting v=2, then we tear it down to sig,cig,hig,norm_base 4 parts. 

$$
\begin{aligned}
cig&=signals(t)*\cos(\phi(t))
sig&=signals(t)*\sin(\phi(t))
hig&=signals(t)
\end{aligned}
$$

norm_base just for normalization

Using suggestion: 
1. filter to 0.1-1Hz before mcc would be nice
2. regarding to feature of CC, total shift should less than 10% of totoal length


## 中文介绍

根据江国明，2012 地球物理学报的内容与脚本改写。

相位加权的多道互相关算法，原式如下 

jupyter notebook里有完整的功能，不需要py文件，但没整理好。

$$ 
C_{\mathrm{mmmcc}}(\tau)=\frac{\sum_{k=1}^{N} S_{j}(t+k \cdot \delta t-\tau) \cdot S_{k}(t+k \cdot \delta t) \cdot\left|\cos \left(\frac{\phi_{j}(t+k \cdot \delta t-\tau)-\phi_{k}(t+k \cdot \delta t)}{2}\right)\right|^{v}}{\sqrt{\sum_{k=1}^{N}\left|S_{j}(t+k \cdot \delta t-\tau)\right|^{2} \cdot \sum_{k=1}^{N}\left|S_{k}(t+k \cdot \delta t)\right|^{2}}}
$$

分为三个部分

相位叠加部分，可以用倍角公式简化，$v=2$时。 

$$ 
\begin{aligned} 
&\cos^2[ \frac{1}{2} \phi_j(t+k\cdot\delta t-\tau) -\frac{1}{2} \phi_j(t+k\cdot\delta t) ]\ 
=& \frac{1}{2} \cos[\phi_j(t-k)]\cos[\phi_k(t)]+ \frac{1}{2} \sin[\phi_j(t-k)]\sin[\phi_k(t)] +\frac{1}{2} 
\end{aligned}
$$

互相关部分

$$
\sum_{k=1}^{N} S_{j}(t+k \cdot \delta t-\tau) \cdot S_{k}(t+k \cdot \delta t)
$$

归一化部分，没什么要解释的。具体的值没有什么意义，要放在整条波形里面。 

$$
\sqrt{\sum_{k=1}^{N}\left|S_{j}(t+k \cdot \delta t-\tau)\right|^{2} \cdot \sum_{k=1}^{N}\left|S_{k}(t+k \cdot \delta t)\right|^{2}}
$$

根据实际经验，建议在匹配前滤波至0.1~1Hz的频段。这时结果比较清晰，高频成分较少。并且实际shift 最好小于波形总长的10%。
