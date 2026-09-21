# Enhanced Tire Brush Model for Vehicle Dynamics Simulation

Summary of Aldo Sorniotti and Mauro Velardocchia, *Enhanced Tire Brush Model for Vehicle Dynamics Simulation*, SAE Technical Paper 2008-01-0595 (2008). This is a summary of the user-supplied PDF, not a transcription of the paper. Page numbers below refer to the 17-page PDF, whose technical paper begins on page 3.

## Purpose and scope

The authors seek a tire model that retains the physical interpretation of a brush model while remaining simple enough for real-time vehicle dynamics calculations. They treat steady-state longitudinal slip, sideslip, and the effects of camber and path curvature. Their main concern is that a model can fit measured total tire forces and aligning moment while describing the contact patch incorrectly. Such a fit would not justify using the model to infer local deflections, pressure distribution, or the effect of changing vertical load. (PDF pp. 3, 12-16.)

## Basic construction

The contact patch is reduced to a line of length $2a$, with coordinate $x$ measured from its center. The authors assume a parabolic normal load per unit length,

$$
q_z(x)=\frac{3F_z}{4a^3}(a^2-x^2), \qquad -a\leq x\leq a,
$$

so integrating $q_z$ over the patch gives the total vertical load $F_z$. A tread brush enters at the leading edge without tangential deflection. While it adheres to the road, longitudinal and lateral deflections grow with the distance traveled through the patch, and local forces follow from separate longitudinal and lateral brush stiffnesses, $c_{pX}$ and $c_{pY}$. Once the local tangential force reaches the available friction, the brush slides. Integrating the local forces and their moment arms gives $F_X$, $F_Y$, and the self-aligning moment $M_Z$. The basic model uses theoretical longitudinal and lateral slip measures derived from slip speed and sideslip angle. (PDF pp. 3-7, especially Eqs. 1-27.)

The extension for camber and curved travel adds a *spin* contribution to lateral brush deflection. The paper treats that contribution as an additional steady-state input, rather than a new differential state. (PDF pp. 7-8, Eqs. 33-38.)

## The enhanced adhesion-to-sliding transition

With unequal longitudinal and lateral brush stiffness, the conventional brush construction changes the *direction* of local force abruptly at the adhesion boundary. Interpreting both sides of that boundary as elastic brush deflections would require an instantaneous jump in deflection. The authors argue that this is not physically plausible. (PDF pp. 5-6 and 8-9, Figs. 3-5 and 7.)

For the usual case $c_{pX}>c_{pY}$, they insert a region between the first loss of adhesion and full two-direction sliding. In that region the brush begins to slide longitudinally, reducing its longitudinal deflection, while its lateral deflection continues to grow as it moves through the patch. The friction limit determines the remaining longitudinal force. A second boundary marks the start of sliding in both directions. The paper also outlines the opposite stiffness ordering, with an approximation for its aligning moment. This change makes the *spatial history of brush deformation through the contact patch* more plausible; it is not a temporal relaxation law for the whole tire. (PDF pp. 8-11, especially Fig. 8 and Eqs. 40-50.)

## Parameters and vertical-load behavior

For practical fitting the authors allow longitudinal and lateral friction to differ, make friction depend on sliding speed and vertical load, reduce lateral brush stiffness with load, and add an approximate carcass-deflection correction to aligning moment. These are modeling and tuning choices, not consequences of the simple brush geometry alone. (PDF pp. 11-12, Eqs. 51-60.)

The paper devotes substantial attention to contact-patch size. A common geometric approximation makes both patch length and width grow approximately with $\sqrt{F_z}$, leaving their ratio nearly fixed. Other hypotheses hold width fixed or infer its load dependence from contact-patch measurements. The authors favor a patch whose length grows faster than its width as load rises. Different choices can give similar forces at one load but substantially different force-versus-load curves. (PDF pp. 13-16, Figs. 14 and 16-21.)

At small slip, the basic model gives approximately

$$
C_{\alpha}=2a^2c_{pY}, \qquad C_{\sigma X}=2a^2c_{pX},
$$

where $a$ is patch half-length and $c_{pX},c_{pY}$ are stiffnesses per unit patch length. Thus changing the patch dimensions changes the predicted small-slip stiffness even if brush properties are held fixed. With their preferred, more realistic load-dependent patch geometry, the authors find that brush stiffness must decrease substantially with increasing vertical load to reproduce the measured less-than-proportional growth of tire stiffness and force. Conversely, a good fit obtained with an unrealistic patch shape can conceal incorrect local mechanics. (PDF pp. 14-16, Eqs. 69-73.)

## Evidence and limits

Figures 10-12 compare predicted longitudinal force, lateral force, and self-aligning moment with experimental tire data at several vertical loads. The paper presents these as examples of steady-state agreement. It does not establish a time-domain relaxation model or validate responses to sudden lift-off and re-contact. The one-dimensional patch, simplified pressure distribution, empirical parameter variations, and approximate carcass correction also limit what can be inferred about actual local tire structure. The authors themselves call for further experimental validation. (PDF pp. 3-4, 12, 16-17.)

## Relevance to Sim3D

This paper is useful when deciding how to make Sim3D's tire force law depend on combined slip and vertical load. In particular, it cautions against calibrating only $F_X$, $F_Y$, and $M_Z$ at one load and then assuming that the contact-patch behavior or load transfer will be correct. Its enhanced brush transition concerns the distribution of forces and deflections *within* a steady contact patch. A separate transient formulation would still be needed if we want tire relaxation states or physically appropriate behavior as the patch shrinks to zero at lift-off. This paragraph is an implication for Sim3D, not a claim that the paper implements those states.
