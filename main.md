---
abstract: |
  Understanding the factors that influence population dynamics across environmental contexts is essential to predict ecological stability. Functional traits influence population growth, which can in turn influence traits and thus create a feedback between population and trait dynamics. We analysed the reciprocal connections between population and trait dynamics in two microbial systems (one autotrophic and one heterotrophic system) and across a range of environmental conditions. When augmenting models of population (trait) dynamics with trait (population density) data, we found no consistent improvement of predictive capacity: this improvement was environment- and system-dependent. Notably, in the cyanobacteria system, models augmented with population density did predict trait dynamics better. In addition, when augmented models were superior, the improvement was substantial. Further investigation of the environmental dependence of trait-growth relationships is recommended.

  **Keywords:** functional traits, density-dependence, trait-dependence, ...
author:
- Mark Holmes, Tessa de Bruin, Nicolas Schtickzelle, Frederik De Laender
bibliography: bibliography.bib
title:
- Environment- and system-specific interactions between traits and population growth
- Supplementary materials
---

# Introduction

Forecasting ecological change requires an understanding of population dynamics, the trajectory of population density through time. A prime predictor of a population's growth rate is its current density. Density-dependence of population growth is generally negative, stimulating growth at low density, and limiting it when density gets high, which effectively puts a cap on population size. Since the classic experiments by Gause (Gause 1934) and Gill (Gill 1972), evidence of negative density-dependence is abundant across biological systems (Sibly et al. 2005). However, means of functional traits across individuals in a population have also been proposed as essential predictors of population growth (Ozgul et al. 2012; Edwards, Litchman, and Klausmeier 2013; Pérez-Ramos et al. 2019; Litchman and Klausmeier 2008; Violle et al. 2007).

One challenge when using traits to predict population dynamics is that traits are not static but dynamic. Rapid trait change is widespread (Ellner, Geber, and Hairston 2011) and, while trait change caused by environmental change is arguably studied more often (Grainger and Levine 2021), traits can also change as a result of population density (Gibert, Grady, and Dell 2022). This can happen when population density is strongly-connected to demographic descriptors (e.g. size structure) (Gillooly et al. 2001; Brown et al. 2004; DeLong et al. 2015; Wieczynski et al. 2021). In addition, the trait value itself can be expected to determine trait change, too. This becomes clear when considering body size, often used as a master trait predicting population dynamics (Litchman and Klausmeier 2008). Energetic and physiological costs associated with body size limit an individual's growth (Frank 2009; Schmidt-Nielsen 1984; Haldane 1926) in much the same way as population size limits population growth. Two candidate predictors of trait change therefore stand out: population density and the trait value itself. Both population growth and trait change therefore represent a dynamical system of coupled variables. This coupling can create feedbacks, and is therefore important to unravel if we are to predict the fate of populations subject to perturbations. One way to formalize this interdependence of population and trait dynamics is as follows (Patel, Cortez, and Schreiber 2018):

$$\begin{split}
        \frac{\text{d}N}{N\text{d}t} & = f(N, T) \\
        \frac{\text{d}T}{\text{d}t} & = g(N, T), \\
    \end{split}
\label{eq:Patel}$$

where $N$ and $T$ are the population density and mean trait value across individuals in a population, and $f$ and $g$ are functions describing the effects of these two predictors on the dynamics of population density and traits. Note that this framework is general; the functions $f$ and $g$ could capture a variety of biological mechanisms, including resource competition, and behavioural, epigenetic, or genetic change. In the context of Eq. [\[eq:Patel\]](#eq:Patel){reference-type="ref" reference="eq:Patel"}, the feedback between population and trait dynamics is captured by the dependence of $f$ on $T$ and of $g$ on $N$. If either dependence is zero, feedback is lost.

Environmental conditions can affect population growth and traits. For example, increases in temperature favour smaller cell size in plankton (Fernández-González and Marañón 2021). Furthermore, interactions between growth and traits can be environment-dependent as well, as environmental factors such as temperature have been shown to modulate how traits affect growth (Brown et al. 2004; Gillooly et al. 2001; Kremer, Thomas, and Litchman 2017). One can therefore expect that feedbacks will be environment-dependent, as well as differing among model systems. In systems where the observed traits have a unequivocal contribution to energy or resource acquisition, these relationships are expected to be stronger. A prime example of such systems are primary producers, where traits making up leaf architecture or pigmentation have a direct link to light and carbon fixation and therefore fitness and population growth (Chin et al. 2023). In such systems, trait information may be essential to correctly predict population dynamics (M. Stomp et al. 2004). How the reciprocal effects (of density on traits) differ between systems is less intuitive, however, and will depend on the capacity and speed of the trait change.

Here, we investigate the reciprocal links between population and trait dynamics in two different microbial systems (ciliates and cyanobacteria) across several environmental conditions. To this end, we performed microcosm experiments where we monitored population densities and trait values over time. Given these data, we first asked whether augmenting a model of density-dependence with traits improved predictions of population growth. Similarly, we then tested whether augmenting a model of trait-dependence with population density improved predictions of trait change. Finally, we asked whether the contributions of traits and density to population and trait dynamics, respectively, depend on the environmental context and/or differ between the two model systems.

We found limited evidence to support general feedbacks between population and trait dynamics. Improvements to predictions of these dynamics by including both trait values and population densities was system- and environment-specific. In cases when the augmented model outperformed the single-state model, the prediction increase was substantial.

# Materials and Methods

## Experiments

### Study systems

We addressed these research questions using a number of unicellular, clonally-reproducing strains from two model systems: seven strains from four ciliate genera and four strains from two ecotypes of the picocyanobacteria genus *Synechococcus* (Farrant et al. 2016; Grébert et al. 2018; Sohm et al. 2016). The ciliate genera represent globally-distributed mixotrophs (Weisse and Montagnes 2022; Lu, Gao, and Weisse 2021) and *Synechoccoccus* makes up a substantial proportion of marine phytoplankton biomass and photosynthesis (Dvořák et al. 2014).

We examined the effects in a range of environmental conditions, represented here by temperature and a pollutant --- both commonly-encountered environmental pressures in natural systems. Specifically, we consider an "altered" environment to have increased temperature, increased concentrations of the herbicide atrazine, or both. Temperature change is a ubiquitously-relevant driver to consider as anthropogenically-driven climate change continues to intensify. Atrazine, present in many ecosystems, represents a broad class of herbicides that act via the inhibition of photosynthesis but which also has negative consequences on the fitness of heterotrophic orgnanisms (Shim et al. 2022). Control conditions are hereafter indicated by "C," temperature changes by "T," atrazine changes by "A," and the combined atrazine and temperature changes by "AT."

We grew the strains in monoculture for sufficient time for them to reach their carrying capacity, enabling us to quantify their intrinsic growth rate and density-dependence. Population densities and trait values were measured at regular, frequent intervals ([\[tab:treatment_table,sec:sampling_times\]](#tab:treatment_table,sec:sampling_times){reference-type="ref" reference="tab:treatment_table,sec:sampling_times"}).

::: {#tab:treatment_table}
| Model system  | Exp. length (days) | No. of strains | Temperature ($\mathrm{^\circ C}$) | Pollution (atrazine, $\mathrm{mg\,L}^{-1}$) |
|:--------------|:-------------------|:---------------|:----------------------------------|:--------------------------------------------|
| Ciliates      | $37$               | $7$            | $\bm{22}, 24$                     | $\bm{0}, 20$                                |
| Cyanobacteria | $10$               | $4$            | $\bm{22}, 24$                     | $\bm{0}, 0.1$                               |

: Experimental design of the two model systems, showing how long the experiment lasted, how frequently they were sampled, the number of strains used, the temperature of the warming treatments in degrees Celsius, and the concentration of the pollutant (atrazine). All combinations of pollution and temperature were included in the experiment design. Treatments in bold indicate the reference conditions.
:::

[\[tab:treatment_table\]]{#tab:treatment_table label="tab:treatment_table"}

### Experimental setup

For the ciliate system, the microcosms consisted of 250 mL square bottles with vented caps (DURAN®) containing 100 mL nutrient medium, inoculated with enough cells to initiate population growth. We obtained the required cell density by a injecting a species-specific volume of 1-week old preculture into new medium (100 L, 5 mL, 10 mL, and 10 mL for *T. thermophila*, *Loxocephalus* sp, *P. caudatum* and *S. teres*, respectively). The microcosms were kept in an incubator with a 10:14h light/dark cycle. The nutrient medium consisted of Chalkey's solution supplemented with crushed protist pellets (Carolina biological supply; 0.55 $\mathrm{gL}^{-1}$), filtered through a 0.22 m filter system (BT50 500 mL) to remove all particulate matter and consecutively bacterized with *Bacillus brevis*, *Brevibacillus subtilis* and *Serratia fonticola*. Following bacterization, the medium was incubated for two days at 20°C on a shaker to allow the bacteria to grow before ciliates were added. During the experiment, we replaced 20% of the medium every third day by replacing 4 mL culture with 4 mL 5x concentrated unbacterized nutrient medium, supplemented with atrazine as needed.

For the cyanobacteria system, the microcosms consisted of 6 mL six-well plates (Sarstedt standard flat 6 well plate, 83.3920) of PCRS-11 Red Sea medium (Rippka et al. 2000), inoculated at 5̃000 cells $\mathrm{mL}^{-1}$ inside temperature-controlled incubators (Lovibond Thermostatic Cabinet). They were cultured in white light (6500K LedAquaristik Sky Bar) in a 12:12h light/dark cycle. The six-well plates were constantly mixed 150 RPM with a VWR mini shaker.

Growth media were adjusted to include the appropriate amount of atrazine and microcosms were kept at a constant temperature throughout. All sampling was undertaken in sterile environments to prevent contamination.

### Sampling procedure

Population densities and traits were measured using procedures specifically designed for each species: imaging for ciliates and flow cytometry for cyanobacteria.

Ciliate microcosms were sampled by digital darkfield imaging using an much updated version of the method described by Pennekamp & Schtickzelle (Pennekamp and Schtickzelle 2013). For each microcosm, at each sample time, a series of pictures was taken of a 810µL sub-sample using a Sony A9 camera equipped with an FE 90 mm F2.8 Macro G OSS lens and a darkfield approach (indirect light and black background). Each series consisted of 10 second burst shots at 10 fps in grey scale (other settings: 1/160s, F11, ISO 6400), resulting in 100 pictures per sample. We then proceeded with automated image analysis using Fiji (Schindelin et al. 2012; Van Rossum 2022) to count, track, and characterize the traits of all moving particles (i.e. cell size, cell shape, movement speed, and movement linearity) in the sample (Pennekamp et al. 2014). The traits were population averages of cell size, cell aspect ratio, linearity of the cell's movement, and movement speed. Cell morphology influences the maintenance requirement and resource uptake rates and the motility traits affect the encounter rate of ciliates for their bacterial resources.

Cyanobacteria microcosms were sampled using a flow cytometer (Guava easyCyte 12HT) that measures cell densities and trait values using several channels, which correspond to specific excitation and fluorescence combinations by lasers and detectors. The flow cytometer channels correspond approximately to cell size and concentration of three photosynthetic pigments present in *Synechococcus*: chlorophyll-a, phycocyanin, and phycoerythrin. At each sampling time, 57 L of distilled water was added to counteract evaporation. Then, 200 L of the 6 mL microcosm was removed and replaced with 200 L of fresh medium. This 200 L sample of the culture was diluted with distilled water to be between 50 and 500 cells $\mathrm{L}^{-1}$, the concentration range for which the flow cytometer gives accurate measurements, and 200 L of the appropriately-diluted sample was put into a 96-well plate and processed by the flow cytometer. The flow cytometer then measures 1000 cells from each well to determine the cell density and measure the trait values. The traits were population averages of several flow cytometer channels, indicative of cell size and concentrations of chlorophyll-a, phycocyanin, & phycoerythrin. Cell size again influences resource uptake and maintenance and the photosynthetic pigments influence the light absorption (Marañón 2015; Maayke Stomp et al. 2007). Sampled debris, dead cells, and doublets were removed from the analysed data using the CyanoFilter and PeacoQC packages in R version 2.3.1 (Olusoji et al. 2021; Emmaneel 2022; R Core Team 2023).

## Analyses

The objectives of our analyses were to investigate whether we could best predict the population per-capita growth rate (hereafter "growth," $\gamma$) from density only, or if augmenting our prediction using trait data improved this prediction in our model systems. Similarly, we tested if augmenting this model using trait data to predict trait change
(hereafter $\tau$) with density data improved predictive capacity.

Since there were multiple traits to incorporate, we first summarised the traits into a single aggregate trait (hereafter "trait") using a principal component analysis, PCA, per strain and treatment. This eliminates the impact of trait covariance and collinearity, and ensures that trait information is represented by a single variable, as is population density. We use only the first component of this PCA to represent the traits, which contains most of the explained variance ([\[fig:cilia_PC_var_explained,fig:cyano_PC_var_explained\]](#fig:cilia_PC_var_explained,fig:cyano_PC_var_explained){reference-type="ref" reference="fig:cilia_PC_var_explained,fig:cyano_PC_var_explained"}).

We calculated, for all combinations of strains and environmental conditions, the growth rate at each time point, $\gamma_t=\Delta^{-1}\textrm{log}(N_{t+\Delta}/N_t)$, where $N$ is population density and $\Delta$ is the time between consecutive sampling times. We then fitted two linear models, for all combinations of strains and environmental conditions, of $\gamma_t$ against $N_t$:

$$\label{eq:pcgr}
\gamma_t = \beta_0 + \beta_1 N_t,$$

$$\label{eq:pcgrExtd}
\gamma_t = \beta_0 + \beta_1 N_t + \beta_2 T_t,$$

where $N_t$ (density, continuous) and $T_t$ (trait value, continuous) are predictors, and the $\beta_i$ are regression coefficients. The first model, Eq. [\[eq:pcgr\]](#eq:pcgr){reference-type="ref" reference="eq:pcgr"}, is the standard model of linear, density-dependent growth. The second model, Eq. [\[eq:pcgrExtd\]](#eq:pcgrExtd){reference-type="ref" reference="eq:pcgrExtd"}, is the model augmented with trait data.

Similarly, for all combinations of strains and environmental conditions, the trait rate of change at each time point $\tau_t$, as $\tau_t=\Delta^{-1}(T_{t+\Delta}-T_{t})$, where $T$ is the trait value, and $\Delta$ is again the time between consecutive sampling times. We again fitted two models, but now using the trait values as the response variables:

$$\label{eq:dT}
\tau_t = \beta_0 + \beta_1 T_t,$$

$$\label{eq:dTExtd}
\tau_t = \beta_0 + \beta_1 T_t + \beta_2 N_t,$$

where $N_t$ (density, continuous) and $T_t$ (trait value, continuous) are predictors, and the $\beta_i$ are regression coefficients. Again, Eq. [\[eq:dTExtd\]](#eq:dTExtd){reference-type="ref" reference="eq:dTExtd"} is the augmented version of Eq. [\[eq:dT\]](#eq:dT){reference-type="ref" reference="eq:dT"}.

With these models at hand, we then tested whether the the augmented model improved predictions of population/trait dynamics by comparing the prediction accuracy of the augmented models relative to the basic models. We compared the AIC of these models to determine which model best predicted the data, transforming these AIC values into Akaike scores representing the probability that the augmented models ([\[eq:pcgrExtd,eq:dTExtd\]](#eq:pcgrExtd,eq:dTExtd){reference-type="ref" reference="eq:pcgrExtd,eq:dTExtd"}) were better than the basic models ([\[eq:pcgr,eq:dT\]](#eq:pcgr,eq:dT){reference-type="ref" reference="eq:pcgr,eq:dT"}) (Wagenmakers and Farrell 2004; Burnham and Anderson 2004). Finally, in conjunction with the likelihood that the augmented model was best, we compared the relative change in predictive accuracy, relative to the single-variable model, using the coefficient of determination (R-squared).

# Results

## Growth and trait change

![Fitted model parameters showing intrinsic population/trait change (left), density-dependence (centre), and trait dependence (right) in seven strains of ciliates. Large rows indicate whether growth or trait change is being examined and subrows indicate whether the coefficients are from the single-variable or augmented model. Points indicate the mean, error bars indicate the standard error, and asterisks indicate parameters are significantly different from zero (vertical black line).](figures/cilia/cilia_td_general.pdf){#fig:cilia_td_general width="0.8\\linewidth"}

![As Fig. [1](#fig:cilia_td_general){reference-type="ref" reference="fig:cilia_td_general"}, but for four strains of the cyanobacteria genus *Synechococcus*.](figures/cyano/cyano_td_general.pdf){#fig:cyano_td_general width="0.8\\linewidth"}

Populations of both systems exhibited mostly nonzero intrinsic growth and negative density-dependence ([\[fig:cilia_td_general,fig:cyano_td_general\]](#fig:cilia_td_general,fig:cyano_td_general){reference-type="ref" reference="fig:cilia_td_general,fig:cyano_td_general"}). Intrinsic growth rates of *Tetrahymena* strains were consistently higher than that of the other strains and *Spirostomum* strains exhibited stronger negative density-dependence (i.e. stronger self-limitation). Environmental conditions had moderate and variable effects on how ciliate populations grew. In the cyanobacteria system, intrinsic growth and density-dependence were comparable across strains ([8](#fig:cyano_pcgr){reference-type="ref" reference="fig:cyano_pcgr"}). However, environmental conditions had greater effects on these parameters; the addition of atrazine substantially reduced intrinsic growth rates (to the point of being negative in cases).

The traits of both systems exhibited generally negative trait-dependence ([\[fig:cilia_td_general,fig:cyano_td_general\]](#fig:cilia_td_general,fig:cyano_td_general){reference-type="ref" reference="fig:cilia_td_general,fig:cyano_td_general"}). In the ciliate system, intrinsic trait change ([9](#fig:cilia_dT){reference-type="ref" reference="fig:cilia_dT"}) was generally not significantly different from zero, but estimates varied across environmental conditions in *Tetrahymena* strains. The trait-dependence of trait change in ciliates was also strongest for these strains. In the cyanobacteria system, intrinsic trait change was often significantly different from zero and mostly positive ([14](#fig:cyano_growth_dtrait){reference-type="ref" reference="fig:cyano_growth_dtrait"}).

## Model selection

![Differences in model performance when adding traits or density to predict growth (top panel) or trait change (bottom panel) different cilia strains (x-axis) in different treatments (colour). The mean values across strains and treatments are shown by the horizontal dashed lines. Values greater than 0.5 indicate that the full model is likely to be better, values less than 0.5 indicate that the single-variable model is likely to better, and values around 0.5 indicate that the two models are approximately equal.](figures/cilia/cilia_AIC.pdf){#fig:cilia_AIC width="\\linewidth"}

In the ciliate system, using the augmented model (including both traits and population densities) to predict growth/trait change generally did not improve predictions ([5](#fig:cilia_AIC_deltaerror){reference-type="ref" reference="fig:cilia_AIC_deltaerror"}). For predictions of growth, only in a few cases (particularly *Tetrahymena* strains and *Loxocephalus* 1) was the augmented model more likely in some environmental conditions. Predictions of trait change with the augmented model were even less likely to outperform the trait-only model.

![As fig. [3](#fig:cilia_AIC){reference-type="ref" reference="fig:cilia_AIC"} but for the cyanobacteria system.](figures/cyano/cyano_AIC.pdf){#fig:cyano_AIC width="0.75\\linewidth"}

In the cyanobacteria system, for treatments without atrazine, the augmented models were generally notably better at predicting the system dynamics ([4](#fig:cyano_AIC){reference-type="ref" reference="fig:cyano_AIC"}). For predictions of population growth, the augmented model was generally likely to best describe the data, with some exceptions in which the added trait data provided no benefit to predictions. For predictions of trait change, the augmented model was more consistently more likely to be better, with only some cases in the treatments with atrazine where it failed to improve the predictive ability.

![Comparing the change in predictive accuracy using the augmented model (both population and traits) as a response to the probability that the augmented model is best for the ciliate system. Model accuracy was calculated using the R-squared statistic. ](figures/cilia/cilia_AIC_deltaerror.pdf){#fig:cilia_AIC_deltaerror width="0.8\\linewidth"}

In the ciliate system, the probability that the augmented model fitted best was more beneficial when making predictions of population growth than trait change ([5](#fig:cilia_AIC_deltaerror){reference-type="ref" reference="fig:cilia_AIC_deltaerror"}). In other words, the cases in which the augmented model was substantially-better fitting (e.g. predictions of growth in *Tetrahymena* 1, control) yielded more accurate predictions of population growth than trait changes.

![As fig. [5](#fig:cilia_AIC_deltaerror){reference-type="ref" reference="fig:cilia_AIC_deltaerror"} but for the cyanobacteria system.](figures/cyano/cyano_AIC_deltaerror.pdf){#fig:cyano_AIC_deltaerror width="0.8\\linewidth"}

In the cyanobacteria system, improvement in predictive accuracy by using the augmented model was substantial only when the augmented model was near-certain to be the best ([6](#fig:cyano_AIC_deltaerror){reference-type="ref" reference="fig:cyano_AIC_deltaerror"}). The augmented model improved predictions of trait change slightly more than predictions of population growth.

We do not observe a consistent, cross-system result as to whether population growth or trait changes are better-predicted by the augmented model. The relationships are clearly environmentally-dependent and system-specific, and variation in augmented model performance is strongly dependent on environmental conditions. That said, for a large number of the cyanobacteria microcosms, particularly in those without atrazine added, the augmented model substantially outperforms the single-variable model ([6](#fig:cyano_AIC_deltaerror){reference-type="ref" reference="fig:cyano_AIC_deltaerror"}).

# Discussion

We have analysed population and trait dynamics in two microbial systems, across various environmental conditions. Our results do not support a general reciprocal relationship between both types of dynamics, but do indicate that trait dynamics are better predicted when considering population density for the cyanobacteria system, in particular. Furthermore, trait effects on growth are dependent on the environment: a given trait value has differently effects on growth, depending on the (a)biotic conditions ([\[fig:cilia_td_general,fig:cyano_td_general\]](#fig:cilia_td_general,fig:cyano_td_general){reference-type="ref" reference="fig:cilia_td_general,fig:cyano_td_general"}) (Gibert, Grady, and Dell 2022; Wieczynski et al. 2021).

Broadly-speaking, the population dynamics of the ciliate system across most environmental conditions were better-described by simple density-dependence, as used in most classical population models e.g. Verhulst (Verhulst 1838). Density-dependent growth is clearly a driving force in these cases as their fitted model coefficients were mostly significantly different from zero ([1](#fig:cilia_td_general){reference-type="ref" reference="fig:cilia_td_general"}) and the augmented model was not strongly supported on the basis of Akaike weights ([3](#fig:cilia_AIC){reference-type="ref" reference="fig:cilia_AIC"}). Trait dynamics were similarly generally better-described without including population density as a covariate ([3](#fig:cilia_AIC){reference-type="ref" reference="fig:cilia_AIC"}). While intrinsic trait dynamics are not as well-founded as the equivalent methods in population dynamics (i.e. there are no real equivalents of trait intrinsic growth or carrying capacity), this basic linear trait-dependent trait change still generally outperformed the augmented model. Trait changes in the ciliate system, then, seem unlikely to depend on the population density and, according to these data, feedbacks between populations and traits in the ciliate system are unlikely.

By contrast, in the cyanobacteria system, the augmented model outperformed the population-only model in predicting population dynamics more often than not and, in the cases where it did not fit well, the population-only model was similarly poorly-fitting ([16](#fig:cyano_R2){reference-type="ref" reference="fig:cyano_R2"}). Furthermore, trait dynamics in the cyanobacteria system were much more consistently well-described by the augmented model. Particularly in the atrazine-free treatments, where the augmented model was near-certain to always be better, but even in the combined atrazine and temperature (AT) treatments, its performance was overall better. Overall, while population is only modestly trait-dependent, the trait dynamics are clearly influenced by population density. In certain cases, therefore, there is potential for feedbacks between populations and traits (e.g. strain V_2524 in the control and temperature treatments).

There are several possible explanations for the patterns observed. Firstly, while we assess the linear relationships between traits and growth rates, it is possible that nonlinear models could more accurately predict the system dynamics. Nonlinear effects of traits on growth are commonplace in natural systems (e.g. the effects of body on growth are generally non-linear (Fernández-González and Marañón 2021)) and such specifics are lost in the simple model formulation here although, for the most part, the linear models used seem appropriate ([\[fig:cilia_pcgr,fig:cilia_dT,fig:cyano_pcgr,fig:cyano_dT\]](#fig:cilia_pcgr,fig:cilia_dT,fig:cyano_pcgr,fig:cyano_dT){reference-type="ref" reference="fig:cilia_pcgr,fig:cilia_dT,fig:cyano_pcgr,fig:cyano_dT"}).

In a similar manner, the better performance of the augmented model in the cyanobacteria system than in the ciliate system may be due to differences in the kind of the traits we measured. The cyanobacteria traits are more directly linked to growth than the ciliate traits. For example, the concentration of photosynthetic pigments within a cell controls absorption of incoming light and therefore growth. Similarly, the effectiveness of light harvesting is likely density-dependent, which makes the pigmentation content subject to plastic change in the time frame of our experiment (M. Stomp et al. 2004). Indeed, high cell densities cause shading, which can alter pigmentation. Pigmentation, therefore may both influence and depend on population growth (Marañón 2015; Maayke Stomp et al. 2007). In contrast, the ciliate traits are as such less-easily linked to growth and thus may require more system-specific analyses. For example, one could summarise the movement and cell size traits into a clearance rate, which may be more strongly-linked to the growth rate of the heterotrophic ciliates (Kiørboe and Saiz 1995; Jacob et al. 2019).

While we have analysed growth curves and trait dynamics to infer coupling between density and traits, akin to time series analysis, a more direct approach would have consisted of a space-for-time approach. In said approach, one experimentally crosses density and trait values, and then monitors short-term growth and trait change (e.g. (Wieczynski et al. 2021)). While manipulating density is feasible, manipulating traits is far less so. One approach could be to grow cultures at different densities, which we know will lead to different traits ([\[fig:cilia_td_general,fig:cyano_td_general\]](#fig:cilia_td_general,fig:cyano_td_general){reference-type="ref" reference="fig:cilia_td_general,fig:cyano_td_general"}), and then create density gradients by dilution. However, such space-for-time approaches do not guarantee a reliable gradient of functional traits and are more experimentally demanding.

We find that the potential for feedbacks is limited in our two experimental systems across the environmental conditions considered. However, in systems with multiple species present (i.e. communities) the potential for feedbacks grows considerably, as indirect interactions may manifest across multiple species' densities and traits (e.g.(Zelnik et al. 2022)). However, one challenge will be that the net result of such feedbacks at the level of population and trait dynamics can be hard to identify, let alone be predicted.

Alternative model-fitting approaches such as autoregressive models may also be used to determine the co-dependence of populations and traits. This may improve the analysis in some ways, e.g. by allowing specific error structures that depend on time. Generally, however, such models require longer time series to be effective and are more of use when forecasting time series than in understanding effects of state variables on each other.

In conclusion, we found that the relationships between populations and functional traits vary with model system and environmental context. While reciprocal interactions between density and traits were apparent for the cyanobacteria system, we did not find broad support for strong feedbacks between population and trait change. More work is needed to unveil the context dependence of growth-trait connections.

# Sampling regime {#sec:sampling_times}

## Ciliate

Due to the large differences in cilate population dynamics, the sampling differed between species. Furthermore, since the most important information required to fit growth curves was towards the start of the experiment (in the exponential growth phase), sampling was more frequent at the beginning of the experiment and then declined in frequency as the populations stabilised. *Tetrahhymena* strains were sampled twice a day during the first week at 09:00 and 17:00. All other species were sampled once a day during the first week at 13:00. During the second week, all species were sample once a day at 13:00. During the third week, all species were sampled every second day at 13:00. Finally, during the fourth and fifth weeks, all species were sampled twice a week on Mondays and Thursdays at 13:00.

## Cyanobacteria

The cyanobacteria microcosms were sampled every day at the same time for each microcosm: between 09:00 and 13:00 depending on the population, as the sampling procedure took some time. Additionally, each microcosm of the cyanobacteria system had three replicates.

# Population density and population change

## Ciliate

![](figures/cilia/cilia_pcgr.pdf){#fig:cilia_pcgr width="0.75\\linewidth"}

## Cyanobacteria

![](figures/cyano/cyano_pcgr.pdf){#fig:cyano_pcgr width="0.75\\linewidth"}

# Trait values and trait change

## Ciliate

![](figures/cilia/cilia_dT.pdf){#fig:cilia_dT width="0.75\\linewidth"}

### Cyanobacteria

![](figures/cyano/cyano_dT.pdf){#fig:cyano_dT width="0.75\\linewidth"}

# Principal component variance explained

## Ciliate

![The explained variance for the first principal component for the cilia models.](figures/cilia/cilia_PC_var_explained.pdf){#fig:cilia_PC_var_explained width="\\linewidth"}

## Cyanobacteria

![The explained variance for the first principal component for the cyanobacteria models.](figures/cyano/cyano_PC_var_explained.pdf){#fig:cyano_PC_var_explained width="\\linewidth"}

# Observed vs predicted trait and population change

## Ciliate

![Observed (x-axis) and predicted (y-axis) population growth (left column) and trait change (left column), including either one (top row) or both predictors (bottom row), for various strains (point types) of ciliates and across environmental conditions (colours).](figures/cilia/cilia_growth_dtrait.pdf){#fig:cilia_growth_dtrait width="\\linewidth"}

## Cyanobacteria

![Observed (x-axis) and predicted (y-axis) population growth (left column) and trait change (left column), including either one (top row) or both predictors (bottom row), for various strains of cyanobacteria (point types) and across environmental conditions (colours).](figures/cyano/cyano_growth_dtrait.pdf){#fig:cyano_growth_dtrait width="\\linewidth"}

# Model goodness-of-fit

## Ciliate

![The coefficient of determination i.e., r-squared, for the cilia system models. The mean values across strains and treatments are shown by the horizontal dashed lines.](figures/cilia/cilia_R2.pdf){#fig:cilia_R2 width="\\linewidth"}

## Cyanobacteria

![As [15](#fig:cilia_R2){reference-type="ref" reference="fig:cilia_R2"} but for the cyanobacteria system models.](figures/cyano/cyano_R2.pdf){#fig:cyano_R2 width="0.833\\linewidth"}

In the ciliate system, the $R^2$ of the model fit was generally below 0.5 for both growth and trait changes for all treatments and was generally unaffected by whether the single- or both-predictor model was used ([15](#fig:cilia_R2){reference-type="ref" reference="fig:cilia_R2"}). *Tetrahymena* strains were generally better-fitting but still without substantial model support. In the cyanobacteria system, the treatments without atrazine were substantially better-fitting than those with atrazine added ([16](#fig:cyano_R2){reference-type="ref" reference="fig:cyano_R2"}). Additionally, while the goodness-of-fit of the growth predictions were only marginally improved by the addition of the full model, the trait change predictions were greatly improved by including the full model.

::: {#refs .references .csl-bib-body .hanging-indent}
::: {#ref-Brown2004 .csl-entry}
Brown, James H., James F. Gillooly, Andrew P. Allen, Van M. Savage, and Geoffrey B. West. 2004. "Toward a Metabolic Theory of Ecology." *Ecology* 85 (7): 1771--89. <https://doi.org/10.1890/03-9000>.
:::

::: {#ref-Burnham2004 .csl-entry}
Burnham, Kenneth P., and David R. Anderson. 2004. "Multimodel Inference: Understanding AIC and BIC in Model Selection." *Sociological Methods & Research* 33 (2): 261--304. <https://doi.org/10.1177/0049124104268644>.
:::

::: {#ref-Chin2023 .csl-entry}
Chin, Alana R. O., Paula Guzmán-Delgado, Anna Görlich, and Janneke HilleRisLambers. 2023. "Towards Multivariate Functional Trait Syndromes: Predicting Foliar Water Uptake in Trees." *Ecology*, no. February: 1--15. <https://doi.org/10.1002/ecy.4112>.
:::

::: {#ref-DeLong2015 .csl-entry}
DeLong, John P., Benjamin Gilbert, Jonathan B. Shurin, Van M. Savage, Brandon T. Barton, Christopher F. Clements, Anthony I. Dell, et al. 2015. "The Body Size Dependence of Trophic Cascades." *The American Naturalist* 185 (3): 354--66. <https://doi.org/10.1086/679735>.
:::

::: {#ref-Dvorak2014 .csl-entry}
Dvořák, Petr, Dale A. Casamatta, Aloisie Poulíčková, Petr Hašler, Vladan Ondřej, and Remo Sanges. 2014. "Synechococcus: 3 Billion Years of Global Dominance." *Molecular Ecology* 23 (22): 5538--51. <https://doi.org/10.1111/mec.12948>.
:::

::: {#ref-Edwards2013 .csl-entry}
Edwards, Kyle F., Elena Litchman, and Christopher A. Klausmeier. 2013. "Functional Traits Explain Phytoplankton Responses to Environmental Gradients Across Lakes of the United States." *Ecology* 94 (7): 1626--35. <https://doi.org/10.1890/12-1459.1>.
:::

::: {#ref-Ellner2011 .csl-entry}
Ellner, Stephen P., Monica A. Geber, and Nelson G. Hairston. 2011. "Does Rapid Evolution Matter? Measuring the Rate of Contemporary Evolution and Its Impacts on Ecological Dynamics." *Ecology Letters* 14 (6): 603--14. <https://doi.org/10.1111/j.1461-0248.2011.01616.x>.
:::

::: {#ref-PeacoQC .csl-entry}
Emmaneel, Annelies. 2022. *PeacoQC: [Peak-based]{.nocase} Selection of High Quality Cytometry Data*. Manual. <https://doi.org/10.18129/B9.bioc.PeacoQC>.
:::

::: {#ref-Farrant2016 .csl-entry}
Farrant, Gregory K., Hugo Doré, Francisco M. Cornejo-Castillo, Frédéric Partensky, Morgane Ratin, Martin Ostrowski, Frances D. Pitt, et al. 2016. "Delineating Ecologically Significant Taxonomic Units from Global Patterns of Marine Picocyanobacteria." *Proceedings of the National Academy of Sciences of the United States of America* 113 (24): E3365--74. <https://doi.org/10.1073/pnas.1524865113>.
:::

::: {#ref-Fernandez-Gonzalez2021 .csl-entry}
Fernández-González, Cristina, and Emilio Marañón. 2021. "Effect of Temperature on the Unimodal Size Scaling of Phytoplankton Growth." *Scientific Reports* 11 (1): 953. <https://doi.org/10.1038/s41598-020-79616-0>.
:::

::: {#ref-Frank2009 .csl-entry}
Frank, S. A. 2009. "The Common Patterns of Nature." *Journal of Evolutionary Biology* 22 (8): 1563--85. <https://doi.org/10.1111/j.1420-9101.2009.01775.x>.
:::

::: {#ref-Gause1934 .csl-entry}
Gause, G. F. 1934. "EXPERIMENTAL ANALYSIS OF VITO VOLTERRA'S MATHEMATICAL THEORY OF THE STRUGGLE FOR EXISTENCE." *Science (New York, N.Y.)* 79 (2036): 16--17. <https://doi.org/10.1126/science.79.2036.16-a>.
:::

::: {#ref-Gibert2022 .csl-entry}
Gibert, Jean P., John M. Grady, and Anthony I. Dell. 2022. "Food Web Consequences of Thermal Asymmetries." *Functional Ecology* 36 (8): 1887--99. <https://doi.org/10.1111/1365-2435.14091>.
:::

::: {#ref-Gill1972 .csl-entry}
Gill, Doublas E. 1972. "Density Dependence and Population Regulation in Laboratory Cultures of Paramecium." *Ecology* 53 (4): 701--8. <https://doi.org/10.2307/1934786>.
:::

::: {#ref-Gillooly2001 .csl-entry}
Gillooly, James F., James H. Brown, Geoffrey B. West, Van M. Savage, and Eric L. Charnov. 2001. "Effects of Size and Temperature on Metabolic Rate." *Science (New York, N.Y.)* 293 (5538): 2248--51. <https://doi.org/10.1126/science.1061967>.
:::

::: {#ref-Grainger2021 .csl-entry}
Grainger, Tess Nahanni, and Jonathan M. Levine. 2021. "Rapid Evolution of Life-History Traits in Response to Warming, Predation and Competition: A Meta-Analysis." *Ecology Letters* n/a (n/a). <https://doi.org/10.1111/ele.13934>.
:::

::: {#ref-Grebert2018 .csl-entry}
Grébert, Théophile, Hugo Doré, Frédéric Partensky, Gregory K. Farrant, Emmanuel S. Boss, Marc Picheral, Lionel Guidi, et al. 2018. "Light Color Acclimation Is a Key Process in the Global Ocean Distribution of *Synechococcus* *Cyanobacteria*." *Proceedings of the National Academy of Sciences* 115 (9): E2010--19. <https://doi.org/10.1073/pnas.1717069115>.
:::

::: {#ref-Haldane1926 .csl-entry}
Haldane, John BS. 1926. "On Being the Right Size." *Harper's Magazine* 152: 424--27.
:::

::: {#ref-Jacob2019 .csl-entry}
Jacob, Staffan, Alexis S. Chaine, Michèle Huet, Jean Clobert, and Delphine Legrand. 2019. "Variability in Dispersal Syndromes Is a Key Driver of Metapopulation Dynamics in Experimental Microcosms." *American Naturalist* 194 (5): 613--26. <https://doi.org/10.1086/705410>.
:::

::: {#ref-Kiorboe1995 .csl-entry}
Kiørboe, Thomas, and Enric Saiz. 1995. "Planktivorous Feeding in Calm and Turbulent Environments, with Emphasis on Copepods." *Marine Ecology Progress Series* 122 (1/3): 135--45. <http://www.jstor.org/stable/24852263>.
:::

::: {#ref-Kremer2017a .csl-entry}
Kremer, Colin T., Mridul K. Thomas, and Elena Litchman. 2017. "Temperature- and Size-Scaling of Phytoplankton Population Growth Rates: Reconciling the Eppley Curve and the Metabolic Theory of Ecology." *Limnology and Oceanography* 62 (4): 1658--70. <https://doi.org/10.1002/lno.10523>.
:::

::: {#ref-Litchman2008 .csl-entry}
Litchman, Elena, and Christopher A. Klausmeier. 2008. "Trait-Based Community Ecology of Phytoplankton." *Annual Review of Ecology, Evolution, and Systematics* 39 (1): 615--39. <https://doi.org/10.1146/annurev.ecolsys.39.110707.173549>.
:::

::: {#ref-Lu2021 .csl-entry}
Lu, Xiaoteng, Yunyi Gao, and Thomas Weisse. 2021. "Functional Ecology of Two Contrasting Freshwater Ciliated Protists in Relation to Temperature." *Journal of Eukaryotic Microbiology* 68 (1): 1--16. <https://doi.org/10.1111/jeu.12823>.
:::

::: {#ref-Maranon2015 .csl-entry}
Marañón, Emilio. 2015. "Cell Size as a Key Determinant of Phytoplankton Metabolism and Community Structure." *Annual Review of Marine Science* 7 (1): 241--64. <https://doi.org/10.1146/annurev-marine-010814-015955>.
:::

::: {#ref-Olusoji2021 .csl-entry}
Olusoji, Oluwafemi D., Jurg W. Spaak, Mark Holmes, Thomas Neyens, Marc Aerts, and Frederik De Laender. 2021. "[cyanoFilter]{.nocase}: An R Package to Identify Phytoplankton Populations from Flow Cytometry Data Using Cell Pigmentation and Granularity." *Ecological Modelling* 460 (November): 109743. <https://doi.org/10.1016/j.ecolmodel.2021.109743>.
:::

::: {#ref-Ozgul2012 .csl-entry}
Ozgul, Arpat, Tim Coulson, Alan Reynolds, Tom C. Cameron, and Tim G. Benton. 2012. "Population Responses to Perturbations: The Importance of Trait-Based Analysis Illustrated Through a Microcosm Experiment." *American Naturalist* 179 (5): 582--94. <https://doi.org/10.1086/664999>.
:::

::: {#ref-Patel2018 .csl-entry}
Patel, Swati, Michael H. Cortez, and Sebastian J. Schreiber. 2018. "Partitioning the Effects of Eco-Evolutionary Feedbacks on Community Stability." *The American Naturalist* 191 (3): 381--94. <https://doi.org/10.1086/695834>.
:::

::: {#ref-Pennekamp2014 .csl-entry}
Pennekamp, Frank, Katherine A. Mitchell, Alexis Chaine, and Nicolas Schtickzelle. 2014. "Dispersal Propensity in Tetrahymena Thermophila Ciliates-a Reaction Norm Perspective." *Evolution; International Journal of Organic Evolution* 68 (8): 2319--30. <https://doi.org/10.1111/evo.12428>.
:::

::: {#ref-Pennekamp2013 .csl-entry}
Pennekamp, Frank, and Nicolas Schtickzelle. 2013. "Implementing Image Analysis in Laboratory-Based Experimental Systems for Ecology and Evolution: A Hands-on Guide." *Methods in Ecology and Evolution* 4 (5): 483--92. <https://doi.org/10.1111/2041-210X.12036>.
:::

::: {#ref-Perez-Ramos2019 .csl-entry}
Pérez-Ramos, Ignacio M., Luis Matías, Lorena Gómez-Aparicio, and Óscar Godoy. 2019. "Functional Traits and Phenotypic Plasticity Modulate Species Coexistence Across Contrasting Climatic Conditions." *Nature Communications* 10 (1): 2555. <https://doi.org/10.1038/s41467-019-10453-0>.
:::

::: {#ref-RCoreTeam2023 .csl-entry}
R Core Team. 2023. *R: A Language and Environment for Statistical Computing*. Manual. Vienna, Austria: R Foundation for Statistical Computing.
:::

::: {#ref-Rippka2000 .csl-entry}
Rippka, R., T. Coursin, W. Hess, C. Lichtle, D. J. Scanlan, K. A. Palinska, I. Iteman, F. Partensky, J. Houmard, and M. Herdman. 2000. "Prochlorococcus Marinus Chisholm Et Al. 1992 Subsp. Pastoris Subsp. Nov. Strain PCC 9511, the First Axenic Chlorophyll A2/B2-Containing Cyanobacterium (Oxyphotobacteria)." *International Journal of Systematic and Evolutionary Microbiology* 50 (5): 1833--47. <https://doi.org/10.1099/00207713-50-5-1833>.
:::

::: {#ref-Schindelin2012 .csl-entry}
Schindelin, Johannes, Ignacio Arganda-Carreras, Erwin Frise, Verena Kaynig, Mark Longair, Tobias Pietzsch, Stephan Preibisch, et al. 2012. "Fiji: An Open-Source Platform for Biological-Image Analysis." *Nature Methods* 9 (7): 676--82. <https://doi.org/10.1038/nmeth.2019>.
:::

::: {#ref-Schmidt-Nielsen1984 .csl-entry}
Schmidt-Nielsen, Knut. 1984. *Scaling: Why Is Animal Size So Important?* Cambridge University Press.
:::

::: {#ref-Shim2022 .csl-entry}
Shim, Kyu Young, Vrinda Sukumaran, In Cheol Yeo, Heesang Shin, and Chang Bum Jeong. 2022. "Effects of Atrazine and Diuron on Life Parameters, Antioxidant Response, and Multixenobiotic Resistance in Non-Targeted Marine Zooplankton." *Comparative Biochemistry and Physiology Part - C: Toxicology and Pharmacology* 258 (May): 109378. <https://doi.org/10.1016/j.cbpc.2022.109378>.
:::

::: {#ref-Sibly2005 .csl-entry}
Sibly, Richard M., Daniel Barker, Michael C. Denham, Jim Hone, and Mark Pagel. 2005. "On the Regulation of Populations of Mammals, Birds, Fish, and Insects." *Science (New York, N.Y.)* 309 (5734): 607--10. <https://doi.org/10.1126/science.1110760>.
:::

::: {#ref-Sohm2016 .csl-entry}
Sohm, Jill A, Nathan A Ahlgren, Zachary J Thomson, Cheryl Williams, James W Moffett, Mak A Saito, Eric A Webb, and Gabrielle Rocap. 2016. "Co-Occurring Synechococcus Ecotypes Occupy Four Major Oceanic Regimes Defined by Temperature, Macronutrients and Iron." *The ISME Journal* 10 (2): 333--45. <https://doi.org/10.1038/ismej.2015.115>.
:::

::: {#ref-Stomp2004 .csl-entry}
Stomp, M., J. Huisman, F. De Jongh, A. J. Veraart, D. Gerla, M. Rijkeboer, B. W. Ibelings, U. I. Wollenzien, and L. J. Stal. 2004. "Adaptive Divergence in Pigment Composition Promotes Phytoplankton Biodiversity." *Nature* 432 (7013): 104--7. <https://doi.org/10.1038/nature03044>.
:::

::: {#ref-Stomp2007a .csl-entry}
Stomp, Maayke, Jef Huisman, Lucas J. Stal, and Hans C. P. Matthijs. 2007. "Colorful Niches of Phototrophic Microorganisms Shaped by Vibrations of the Water Molecule." *ISME Journal* 1 (4): 271--82. <https://doi.org/10.1038/ismej.2007.59>.
:::

::: {#ref-VanRossum2022 .csl-entry}
Van Rossum, G. 2022. "The Python Language Reference: Expressions." *Python Reference Manual*, 11.
:::

::: {#ref-Verhulst1838 .csl-entry}
Verhulst, Pierre-François. 1838. "Notice Sur La Loi Que La Population Suit Dans Son Accroissement." *Correspondence Mathematique Et Physique* 10: 113--29.
:::

::: {#ref-Violle2007a .csl-entry}
Violle, Cyrille, Marie-Laure Navas, Denis Vile, Elena Kazakou, Claire Fortunel, Irène Hummel, and Eric Garnier. 2007. "Let the Concept of Trait Be Functional!" *Oikos* 116 (5): 882--92. <https://doi.org/10.1111/j.0030-1299.2007.15559.x>.
:::

::: {#ref-Wagenmakers2004 .csl-entry}
Wagenmakers, Eric-Jan, and Simon Farrell. 2004. "AIC Model Selection Using Akaike Weights." *Psychonomic Bulletin & Review* 11 (1): 192--96. <https://doi.org/10.3758/BF03206482>.
:::

::: {#ref-Weisse2022 .csl-entry}
Weisse, Thomas, and David J. S. Montagnes. 2022. "Ecology of Planktonic Ciliates in a Changing World: Concepts, Methods, and Challenges." *Journal of Eukaryotic Microbiology* 69 (5): 1--16. <https://doi.org/10.1111/jeu.12879>.
:::

::: {#ref-Wieczynski2021 .csl-entry}
Wieczynski, Daniel J, Pranav Singla, Adrian Doan, Alexandra Singleton, Ze-yi Han, and Samantha Votzke. 2021. "Linking Species Traits and Demography to Explain Complex Temperature Responses Across Levels of Organization." *Proceedings of the National Academy of Sciences of the United States of America* 118 (42). <https://doi.org/10.1073/pnas.2104863118/-/DCSupplemental.Published>.
:::

::: {#ref-Zelnik2022a .csl-entry}
Zelnik, Yuval R., Nuria Galiana, Matthieu Barbier, Michel Loreau, Eric Galbraith, and Jean-François Arnoldi. 2022. "How Collectively Integrated Are Ecological Communities?" Preprint. Ecology. <https://doi.org/10.1101/2022.12.29.522189>.
:::
:::
