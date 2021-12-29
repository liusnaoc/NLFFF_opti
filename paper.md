---
title: 'NLFFF: A Python package for a 3D non-linear force-free magnetic field extrapolation in cartesian coordinate system'
tags:
  - Python
  - astronomy
  - solar physics
  - magnetic field
  - force-free field
authors:
  - name:  Liu S  # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-0872-7098
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: National Astronomical Observatories, CAS, Beijing, 100101, China
   index: 1
 - name: School of Astronomy and Space Sciences, University of CAS, Beijing, 100101, China
   index: 2

date: 29 December 2021
bibliography: paper.bib

---

# Summary
At present, the measurement of the solar magnetic field is limited to the photosphere and chromosphere.
Alternatively, coronal magnetic field is often obtained by magnetic field extrapolation under 
some theoretical model. Force-free model is a common accepted one for corona magnetic field 
extrapolation due the low $\beta$ ($\beta$ is the ratio of gas to magnetic press) plasma environment. 
Force-free magnetic field follows:  $\nabla\times\textbf{B}=\alpha\textbf{B}$ and $\nabla\cdot\textbf{B}=0$.
Because Python has become the most popular and open source programming language within
the solar physics and astronomy community (Bobra et al., 2020).
Hence,a python package performs a non-linear force-free-field (NLFFF) extrapolation of the solar magnetic field using the
optimization method of of d of of Wheatland et al. (2000) is developed. The method optimize an objective function, L,
which is the integral of $|\textbf{J}\times\textbf{B}|^2+|\nabla\cdot\textbf{B}|^2$ ($\textbf{J}=(\nabla\times\textbf{B})/\mu$) over a volume to be minimum. The programe can extrapolated 3D corona magnetic field in cartesian coordinate system satisfying force-free model from the observation of photosphere magnetic field as boundary. Here the force-free field can be the potential, linear and non-linear force-free field.
The codes are referred to and edited from SSW package, which can give preliminary NLFFF. For more accurate results, one can alternatively use special calculation program with high speed and high precision results.

# Statement of need
The magnetic field in the solar active active region dominate most solar activities such as filaments eruptions, flares and coronal mass ejections (CMEs). All these eruptions phenomena are energetic events due to explosive release of magnetic field energy. Hence,
the knowledge of magnetic field is necessary to understand solar activities. Magnetic fields fill in the whole spaces above the active regions, their evolutions should be exhibited in whole spaces of active regions. However the magnetic fields with sufficient resolution and accuracy limited to be measured on the photosphere, while magnetic field measurements in chromosphere and corona are only available for a few special cases. To study the properties of spatial magnetic fields, the magnetic field extrapolations usually be employed in practice. Force free extrapolation is regarded as an classical method used to study solar spatial magnetic field
recently. Force free extrapolations are based on the assumption that there is no Lorentz force in the whole space of active region,
which can be expressed by $(\nabla\times\textbf{B}\times\textbf{B} = 0$ mathematically. The spatial/coronal magnetic fields can be reconstructed from this physical assumption (namely, $(\nabla\times\textbf{B}=\alpha\textbf{B}$ and $\nabla\cdot\textbf{B} = 0$), in which the observed photospheric magnetic fields are taken as a boundary conditions. There are lots of methods can be applied to extrapolate force-free field, for potential field $\alpha=0$ and linear force-free $\alpha=constant$ the extrapolation methods are well developed. For non-linear force-free extrapolation, the methods are developing with varioius algorithms under the process of development. Here, the program using Python is developed, which can obtain the potential and linear force-free field using lassical algorithms and non-linear force-free field using optimization algorithm.

# Mathematics
Green's function calculation of Chiu \& Hilton (1977) are used for potential 
and linear force-free field extrapolation.

\begin{equation}
\begin{split}
B_{i}= \frac{1}{2\pi}  \int_{0}^{L_{x}}
\int_{0}^{L_{y}}dx^{'}dy^{'}G_{i}(x,y,z;x^{'},y^{'})B_{z}(x^{'},y^{'},0),
\hspace{2ex} i=x,y,z,
\end{split}
\end{equation}

\begin{equation}
\begin{split}
G_{x}=\frac{x-x^{'}}{R}\frac{\partial\Gamma}{\partial z}+\alpha\Gamma \frac{y-y^{'}}{R},\\
G_{y}=\frac{y-y^{'}}{R}\frac{\partial\Gamma}{\partial z}-\alpha\Gamma \frac{x-x^{'}}{R},\\
G_{z}=-\frac{\partial\Gamma}{\partial R}-\frac{\Gamma}{R},\\
\end{split}
\end{equation}
\begin{equation} \Gamma=\frac{z}{Rr}cos(\alpha
r)-\frac{1}{R}cos(\alpha z),
\end{equation}
\begin{equation}
R=\sqrt{(x-x^{'})^{2}+(y-y^{'})^{2}}.
\end{equation}
Where $B_{z}(x^{'},y^{'},0)$ is the light of sight magnetic field of photosphere. 
Figures \autoref{fig:potent}  and \autoref{fig:linear}  show the field lines distribution of potential and linear force-free extrapolated 3D magnetic field using one of Low and Lou (1990) analytical solution as boundary condition.

Optimization method of Wheatland et al. (2000) are used for non-linear force-free field extrapolation.

 \begin{equation}
 \label{opt1}
 L=\int_{V}[B^{-2}|(\nabla \times \mathbf{B})\times \mathbf{B}|^{2}+|\nabla\cdot\mathbf{B}^{2}|]dV,
\end{equation}

 \begin{equation}
 \label{opt1}
 \dfrac{dL}{dt}=-2\int \frac{\partial \mathbf{B}}{\partial t}\cdot\mathbf{F}dV-\int \frac{\partial \mathbf{B}}{\partial t}\cdot\mathbf{G}d\mathbf{S}
\end{equation}

 \begin{equation}
 \label{opt1}
F=\nabla\times(\mathbf{\Omega}\times\mathbf{B})-\mathbf{\Omega}\times(\nabla\times\mathbf{B})-\nabla(\mathbf{\Omega}\cdot\mathbf{B})+\mathbf{\Omega}(\nabla\cdot\mathbf{B})+\mathbf{\Omega}^{2}\mathbf{B}
\end{equation}

 \begin{equation}
 \label{opt1}
G=\mathbf{n}\times(\mathbf{\Omega}\times\mathbf{B})-\mathbf{n}(\mathbf{\Omega}\cdot\mathbf{B})
\end{equation}


 \begin{equation}
 \label{opt1}
\mathbf{\Omega}=B^{-2}[(\nabla\times\mathbf{B})\times\mathbf{B}-(\nabla\cdot\mathbf{B})\mathbf{B}]
\end{equation}


\begin{equation}\label{lbf}
 \frac{\partial \mathbf{B}}{\partial t}=\mu F
\end{equation}


\begin{equation}\label{lbf}
 \frac{\partial \mathbf{B}}{\partial t}=0
\end{equation}

\begin{equation}\label{lbf}
 \frac{dL}{dt}=-2\int_{V}\mu F^{2}dV.
\end{equation}
Where $t$ is a introducing artificial pseudo-time parameter, for the optimization iterative algorithm.
Figures \autoref{fig:nonlinear} shows the field lines distribution of extrapolated 3D magnetic field using one of Low \& Lou (1990) analytical solution as boundary condition.



# Figures

[Field lines distribution of potential magnetic field extrapolated, $\alpha=0$.\label{fig:potent}](potent.eps)

[Field lines distribution of linear force-free magnetic field extrapolated, $\alpha=<\alpha_{boundary}>$.\label{fig:linear}](linear.eps)

[Field lines distribution of non-linear force-free magnetic field extrapolated.\label{fig:nonlinear}](nonlinear.eps)




# Acknowledgements
Support from scientific research projects of (Grant No. XDA15320301, XDA15320302, XDA15052200),
 and acknowledge contributions of NLFFF in SSW package.


\begin{thebibliography}{}


Chiu, Y.T. \& Hilton, H.H., 1977, ApJ., \textbf{212},873-885.
\href{http://dx.doi.org/10.1086/155111}{DOI}.
\href{https://ui.adsabs.harvard.edu/abs/1977ApJ...212..873C}{ADS}.\\


Bobra, M.G., Mumford, S.J., Hewett, R.J., Christe, S.D., Reardon, K., Savage, S., Ireland, J., Pereira, Tiago M.D., Chen, B., Prez-Suarez, D., 2020, SolPhys., \textbf{295},57.
\href{http://dx.doi.org/10.1007/s11207-020-01622-2}{DOI}.
\href{https://ui.adsabs.harvard.edu/abs/2020SoPh..295...57B/}{ADS}.\\




Wheatlan, M.S., Sturrock, P.A. \& Roumeliotis, G.,
2000, ApJ., \textbf{540}, 1150-1155.
\href{http://dx.doi.org/10.1086/309355}{DOI}.
\href{https://ui.adsabs.harvard.edu/abs/2000ApJ...540.1150W}{ADS}.



\end{thebibliography}

