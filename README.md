# QTAG Basis with Semi-local GBF Momentum

As in all GBF methods, we start with a wavefunction represented as a superposition of $$\( N_b \)$$ Gaussian functions, $$\( \{g_n\} \)$$:

$$g_n(x, t) := \left(\frac{a_n}{\pi}\right)^{1/4} \exp\left(-\frac{a_n (x - q_n)^2}{2} + i p_n (x - q_n) + i s_n\right)$$

The wavefunction can be expressed as:

$$\psi(x, t) = \sum_{n=1}^{N_b} c_n(t) g_n(x; \lambda_n(t))$$

where $$\( \{c_n\} \)$$ are the expansion coefficients.  $$\( \lambda_n(t) \)$$ is a vector of the time-dependent parameters of the $$\( n \)-th$$ GBF, $$\( g_n(x, t) \)$$, which in QTAG includes the center position $$\( q_n \)$$, the momentum $$\( p_n \)$$ defining the linear phase of $$\( g_n \)$$, and the real width parameter $$\( a_n \)$$. Therefore, 

$$\lambda_n(t) = (q_n, p_n, a_n)$$


The EOMs for the QTAG parameters are:

$$p_n := \Im\left(\frac{\nabla \psi}{\psi}\right) \bigg|_{x=q_n} $$

$$\frac{dq_n}{dt} = \frac{p_n}{m}$$

$$\frac{da_n}{dt} = -\frac{\nabla p_n}{m} a_n$$

The EOMs for the expansion coefficients, $$\( \{c_n\} \)$$, are:

$$i S \dot{c} = \left(H - i \dot{S}\right)c$$

where $$\( H \)$$ and $$\( S \)$$ are the Hamiltonian and the basis overlap matrices, respectively:

$$H_{kn} = \langle g_k | \hat{H} | g_n \rangle$$

$$S_{kn} = \langle g_k | g_n \rangle$$

The non-Hermitian matrix $$\( \dot{S} \)$$ contains the full time-derivative of the $$\( n \)-th$$ basis function with respect to its parameters enumerated with the subscript $$\( l \)$$:

$$\dot{S}_{kn} = \langle g_k | \dot{g_n} \rangle$$

$$\dot{g_n} = \sum_{l} \frac{\partial g_n}{\partial \lambda_{n,l}} \dot{\lambda}_{n,l}$$

In QTAG, we do not solve the Bohmian EOM for the momentum of a QT at \( q \):

$$\frac{d}{dt} p_{\text{bohm}} := -\nabla(V + U) \bigg|_{x=q}$$

which involves a possibly singular gradient of the quantum potential $$\( U(x, t) \)$$:

$$U(x, t) := -\frac{1}{2m} \frac{\nabla^2 |\psi|}{|\psi|} $$

formally associated with the division by $$\( |\psi| \)$$. Instead, the GBF parameters are updated according to $$\( p_n \)$$ and its gradient, defined by $$\( \psi(x, t) \)$$ expressed in a basis.

Even though $$\( \psi(x, t) \$$ is available analytically as a basis expansion, the straightforward implementation may generate large values near the wavefunction nodes, leading to numerically unstable trajectory dynamics.

Below we outline a new procedure for computing $$\( p_n \)$$, in which division by $$\( \psi \)$$ is circumvented and the "spatial resolution" of $$\( p_n \)$$ is controlled by a single parameter, $$\( \beta \)$$, defined below, controlling the locality of the fitting.

The action of the momentum operator on $$\( \psi(x, t) \)$$ can be expressed as:

$$\hat{p} \psi(x, t) = -i \eta(x, t) \psi(x, t)$$



where $$\( \eta(x, t) \)$$ is a complex function:

$$\eta(x, t) := \frac{\nabla \psi(x, t)}{\psi(x, t)}$$

whose real and imaginary components, $$\( r(x, t) \)$$ and $$\( p(x, t) \)$$, are referred to as the non-classical $$\( r(x, t) \)$$ and classical $$\( p(x, t) \)$$ momentum components:

$$r := \frac{\nabla |\psi|}{|\psi|}, \quad p := \nabla (\arg \psi) $$

Let us focus on fitting the complex function $$\( \eta \)$$, rather than $$\( p \)$$. We will define the GBF momentum $$\( p_n \)$$ (and if desired $$\( r_n \))$$ from the minimization of the $$\( q_n \)$$-centered error functional:

$$I_n = \int \frac{|\eta - \tilde{\eta}_n|^2}{|\psi|^2} w_n dx$$

where $$\( w_n \)$$ is a Gaussian (or other spatially localized) window function:

$$w_n(x; q_n, \beta) := \sqrt{\frac{\beta}{\pi}} \exp\left(-\beta (x - q_n)^2\right) $$

In particular, the linear fit centered at $$\( q_n \)$$, \( \tilde{\eta}_n = d_1 + d_2(x - q_n) \), directly yields $$\( p_n \)$$ and $$\( \nabla p_n \)$$ – the latter needed to update the GBF width – as the imaginary parts of the optimal $$\( d_1 \)$$ and $$\( d_2 \)$$, respectively. These fitting parameter values, arranged as a vector:

$$d = \begin{pmatrix} d_1 \\ d_2 \end{pmatrix}$$

solve the following system of linear equations for each $$\( q_n \)$$:

$$M d = -b $$

$$b_j = \langle w_n \psi | (x - q_n)^{j-1} | \nabla \psi \rangle$$

The following equations define the system of linear equations for computing the fitting parameters:

The equation for $$\( M \)$$ and $$\( b \)$$:
$$M d = b \$$

The equation for $$\( b_j \)$$:
$$b_j = \langle w_n \psi | (x - q_n)^{j-1} | \nabla \psi \rangle $$

The equation for \( M_{ij} \):
$$M_{ij} = \langle w_n \psi | (x - q_n)^{i+j-2} | \psi \rangle, \quad i, j \in \{1, 2\}.$$

Using $$\( |\psi|^2 \)$$ as the weighting function in $$\( I_n \)$$ removes the division by $$\( \psi \)$$ in the integrand. Additionally, the localization of $$\( w_n \)$$ tunes the resolution of $$\( \tilde{\eta}_n \)$$ from low $$(\( \beta \to 0 \))$$ to high $$(\( \beta \to \infty \))$$, transitioning from a global to a local fit.

Reference: <span style="color:blue;">Garashchuk, S. and Großmann, F., 2024.</span> <em>Assessing the Accuracy of Quantum Dynamics Performed in the Time-Dependent Basis Representation.</em> <span style="color:blue;">The Journal of Physical Chemistry A.</span>



