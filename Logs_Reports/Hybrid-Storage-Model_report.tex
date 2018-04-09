\documentclass{article}
\usepackage[utf8]{inputenc}

\title{Hybrid Storage Model}
\author{Alok Deshpande}
\date{March 11, 2018}

%\usepackage{natbib}
\usepackage{graphicx}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{mathtools}
\usepackage{bm}

%\usepackage[backend=biber]{biblatex}
%\addbibresource{references.bib}

\DeclareMathOperator{\E}{\mathbb{E}}

\begin{document}
	
	\maketitle
	
	\section{Introduction}
	Battery life is a factor that is critical to the operation of electric vehicles. Amongst its many effects on vehicle performance, it dictates the maximum distance that can be travelled on a single charging. Naturally, one would hope to maximize battery life.
	
	In doing so, a common strategy being applied by vehicle designers is using a hybrid power supply, composed of a combination of a battery and a secondary energy storage. \cite{thounthong2009energy} The reason is that each of these storage devices is best suited for a specific energy demand (load) profile for the vehicle's operations, which may vary during the course of the trip. By using a combination of storages, one may improve the efficiency of the energy delivery, and hence the battery life. Note that the secondary storage is usually a supercapacitor and so will be referred to as such in this paper, even though the results would apply for other devices too.
	
	In order to improve its lifetime, experiments show that the battery should be operated with constant, low output power. That being said, it has greater total energy capacity. The supercapacitor, on the other hand, has the opposite characteristics: high power output possible, but low energy capacity. Hence, it is clear that optimally discharging energy from the battery and the supercapacitor is important in order to maximize the battery lifetime within the constraint of always supplying the demand, which varies over time depending on vehicle operation.
	
	---
	
	This model describes the energy flows in a two-storage system. It is formulated as a dynamic programming problem.
	
	\section{Model}
	
	\subsection{Definitions}
	\begin{itemize}
		\item Indices\\
		$t$: discrete time step index\\
		$j$: storage device number (1: battery, 2: supercapacitor)\\
		\item Parameters\\
		$\alpha^{C}$: charging efficiency\\
		$\alpha^{D}$: discharging efficiency\\
		$\beta$: storage efficiency factor (constant)\\
		$N$: number of steps (DP horizon)\\
		$K$: cost weighting factor for rate (relative to cost of energy loss)\\
		\item Variables\\
		$L$: load energy demand (random variable)\\
		$E$: energy state of storage device\\
		$D$: energy released by discharging (AFTER loss)\\
		$C$: energy consumed by charging (BEFORE loss)\\
		$J$: cost function
	\end{itemize}
	
	NOTE: $C_{1}$ does not exist because not possible to charge the battery while driving. (Assuming no regenerative braking at the moment.)
	
	\subsection{Constraints}
	\begin{itemize}
		\item Supply-demand balance: 
		\begin{equation} \label{eq:BalanceEqn}\left[D_{1}(t)\right] + \left[D_{2}(t) - C_{2}(t)\right] = L(t) \end{equation}
		
		\item Bounds on stored energy: 
		\begin{equation}E_{j}^{min}\leq E_{j}(t)\leq E_{j}^{max}\end{equation}
		\item Bounds on charging:
		\begin{equation}0\leq C_{2}(t)\leq C_{2}^{max}\end{equation}
		\item Bounds on discharging:
		\begin{equation}0\leq D_{j}(t)\leq D_{j}^{max}\end{equation}
	\end{itemize}
	
	\subsection{Recursive State Equations}
	The state of the system is the energy in a storage device ($E_{j}(t)$). This evolves according to the charging and discharging of the storage device, which is the control.
	
	The following recursive equations describe the changes in the state, including due to constant leakage loss:
	
	\begin{equation}\label{eq:StateEq1}E_{1}(t+1)=\beta_{1}E_{1}(t)+\left[-\frac{1}{\alpha_{1}^{D}}D_{1}(t)\right] \end{equation}
	
	\begin{equation}\label{eq:StateEq2}E_{2}(t+1)=\beta_{2}E_{2}(t)+\left[\alpha_{2}^{C}C_{2}(t)-\frac{1}{\alpha_{2}^{D}}D_{2}(t)\right]\end{equation}
	
	Substituting \eqref{eq:BalanceEqn} into \eqref{eq:StateEq2}, one obtains:
	
	\begin{equation}\label{eq:DPStateEq2}E_{2}(t+1)=\beta_{2}E_{2}(t)+\left[\alpha_{2}^{C}[D_{1}(t)+D_{2}(t)-L(t)]-\frac{1}{\alpha_{2}^{D}}D_{2}(t)\right]\end{equation}
	
	
	\subsection{Cost Function}
	\begin{itemize}
		\item Minimize discharge rate
		for the first storage device (battery):
		\begin{equation}J_{rate}=min\left[\sum_{t=0}^{N-1}K\left[D_{1}(t)\right]^{2}\right]\end{equation}
		This is convex.
		\item Minimize power loss due to energy transfers
		\begin{equation}J_{loss}=min\left[\sum_{t=0}^{N-1}
		(1-\alpha_{1}^{D})D_{1}(t)+
		(1-\alpha_{2}^{C})C_{2}(t)+
		(1-\alpha_{2}^{D})D_{2}(t)
		\right]\end{equation}
		Substituting \eqref{eq:BalanceEqn}, one obtains:
		\begin{equation}J_{loss}=\mathop{\E}_{L(t)} \Biggl\{min\left[\sum_{t=0}^{N-1}
		(1-\alpha_{1}^{D})D_{1}(t)+
		(1-\alpha_{2}^{C})[D_{1}(t)+D_{2}(t)-L(t)]+
		(1-\alpha_{2}^{D})D_{2}(t)
		\right]\Biggr\}\end{equation}
		This is constrained to be non-negative.
	\end{itemize}
	
	Hence the combined cost function is convex:
	\begin{equation}J=\mathop{\E}_{\substack{ L(t) \\ \forall t\in{1\dots N} }} \Biggl\{min\left[\sum_{t=0}^{N-1}K\left[D_{1}(t)\right]^{2} + (1-\alpha_{1}^{D})D_{1}(t)+
	(1-\alpha_{2}^{C})[D_{1}(t)+D_{2}(t)-L(t)]+
	(1-\alpha_{2}^{D})D_{2}(t)\right]\Biggr\}\end{equation}
	It is chosen to take the expectation after the minimization. This is done so that the net energy discharged (control, $D$) exactly matches the demand, $L$, at all times, and not its expected value.
	
	In the general case, for state $x(t)$, control $u(t)$ and random perturbation $w(t)$, the cost function may be expressed as:
	
	\begin{equation}J=\mathop{\E}_{w(t)} \Biggl\{\min_{u}\left[\sum_{t=0}^{N-1}g(x(t),u(t),w(t))\right]\Biggr\}\end{equation}
	where $g(\cdot)$ is the stage cost.
	
	This can be re-written in the form of Bellman's equation, which allows the problem to be solved by recursion:
	\begin{multline}
	J_{t}[x(t),w(t)]=\min_{u} g(x(t),u(t),w(t)) + \mathop{\E}_{w(t+1)} \{J_{t+1}[f(x(t),u(t),w(t)),w(t+1)]\}
	\end{multline}
	Note that by choosing the perturbation $w(t)$ to be a component of the state as well, can determine the optimal control for any arbitrary load at time $t$, in addition to arbitrary state.
	
	Based on this formulation, the optimization problems for the battery and supercapacitor are individually as follow:
	
	\begin{itemize}
		\item Battery storage:\\
		\begin{multline}
		J_{t}[E_{1}(t),L(t)] = \min_{D_{1},D_{2}}
		(1-\alpha_{1}^{D})D_{1}(t) 
		+ K[D_{1}(t)]^{2}
		+(1-\alpha_{2}^{D})[D_{2}(t)]\\	  +(1-\alpha_{2}^{C})[D_{1}(t)+D_{2}(t)-L(t)]
		+\mathop{\E}_{L(t+1)}\{J_{t+1}[f_{1}(E_{1}(t),D_{1}(t)),L(t+1)]\}
		\end{multline}
		
		where $f_{1}(\cdot)$ is \eqref{eq:StateEq1}, the state equation for the battery.
		
		\item Supercapacitor storage:\\
		\begin{multline}
		J_{t}[E_{2}(t),L(t)] = \min_{D_{1},D_{2}}
		(1-\alpha_{1}^{D})D_{1}(t) 
		+ K[D_{1}(t)]^{2}
		+(1-\alpha_{2}^{D})[D_{2}(t)]\\	  +(1-\alpha_{2}^{C})[D_{1}(t)+D_{2}(t)-L(t)]
		+\mathop{\E}_{L(t+1)} \{J_{t+1}[f_{2}(E_{2}(t),D_{1}(t),D_{2}(t),L(t)),L(t+1)]\}
		\end{multline}
		
		where $f_{2}(\cdot)$ is \eqref{eq:DPStateEq2}, the state equation for the supercapacitor.
		
	\end{itemize}
	
	Combining the above gives the final form of the optimization problem of interest:
	\begin{multline}
	J_{t}[E_{1}(t),E_{2}(t),L(t)] = \min_{D_{1},D_{2}}
	(1-\alpha_{1}^{D})D_{1}(t) 
	+ K[D_{1}(t)]^{2}\\
	+(1-\alpha_{2}^{D})[D_{2}(t)]	  +(1-\alpha_{2}^{C})[D_{1}(t)+D_{2}(t)-L(t)]\\
	+\mathop{\E}_{L(t+1)}\{J_{t+1}[f_{1}(E_{1}(t),D_{1}(t)), f_{2}(E_{2}(t),D_{1}(t),D_{2}(t),L(t)), L(t+1)]\}
	\end{multline}
	
	
	
	\section{Finite Horizon DP}
	\begin{multline} \label{eq:FHDP_eq}
	J_{t}[E_{1}(t),E_{2}(t),L(t)] = \min_{D_{1},D_{2}}
	(1-\alpha_{1}^{D})D_{1}(t) 
	+ K[D_{1}(t)]^{2}\\
	+(1-\alpha_{2}^{D})[D_{2}(t)]	  +(1-\alpha_{2}^{C})[D_{1}(t)+D_{2}(t)-L(t)]\\
	+\mathop{\E}_{L(t+1)}\{J_{t+1}[f_{1}(E_{1}(t),D_{1}(t)), f_{2}(E_{2}(t),D_{1}(t),D_{2}(t),L(t)), L(t+1)]\}
	\end{multline}
	\\
	s.t. 
	\begin{itemize}
		\item Bounds on stored energy: 
		\begin{math}E_{j}^{min}\leq E_{j}(t)\leq E_{j}^{max}\end{math}
		\item Bounds on charging:
		\begin{math}0\leq C_{j}(t)\leq C_{j}^{max}\end{math}
		\item Bounds on discharging:
		\begin{math}0\leq D_{j}(t)\leq D_{j}^{max}\end{math}
	\end{itemize}
	where $f_{1}$ and $f_{2}$ are the recursive state equations:\\
	\begin{itemize}
		\item \begin{math}f_{1}(E_{1}(t),D_{1}(t))=\beta_{1}E_{1}(t)+\left[-\frac{1}{\alpha_{1}^{D}}D_{1}(t)\right]\end{math}\\
		\item \begin{math}f_{2}(E_{2}(t),D_{1}(t),D_{2}(t),L(t))=\beta_{2}E_{2}(t)+\left[\alpha_{2}^{C}[D_{1}(t)+D_{2}(t)-L(t)]-\frac{1}{\alpha_{2}^{D}}D_{2}(t)\right]\end{math}
	\end{itemize}
	
	%	\begin{math} \min_{D_{1},D_{2}} [D_{1}+D_{2}-L] \end{math} s.t. \begin{math} constraints \end{math}
	
	\begin{math} D_{1}(t) + D_{2}(t) = L(t) + C_{2}(t) \geq L(t) \end{math}
	
	
	
	\section{Infinite Horizon DP}
	\begin{multline}
	J^{*}[x(N-2),w(N-2)]=\min_{u} g(D_{1}(N-2),D_{1}(N-2),L(N-2))\\
	+\mathop{\E}_{w(N-1)} \{\alpha J^{*}[f_{1}(E_{1}(N-2),D_{1}(N-2)), f_{2}(E_{2}(N-2),D_{1}(N-2),D_{2}(N-2),L(N-2)),L(N-1)]
	\}
	\end{multline}
	subject to the same constraints as above.\\\\
	Here $J^{*}$ is the converged cost function and $0<\alpha<1$ is the discount factor. In the limit, terms with coefficients $\alpha^{N},N\to\infty$ become negligible.
	
	
	
	
	\section{Linear Programming formulation of IHDP}
	Let:
	\begin{itemize}
		\item $X$ be the set of all possible discrete states, indexed by $1\leq i \leq N$
		\item $\lambda(x_{i})$ be the cost of state $x_{i}$
		\item $U$ be the set of all possible discrete controls, indexed by $1\leq p \leq P$
		\item $p_{i,j}(u)$ be the probability of transitioning from state $x_{i}$ to $x_{j}$ by control $u$
		\item $g(x_{i},u_{p})$ be the stage cost
	\end{itemize}
	
	By applying the Bellman operator with discount factor, $T$, repeatedly, the cost will converge to a maximum of $\boldsymbol{J^{*}(x)}$. Hence, the following is a general LP formulation of IHDP, in standard LP form:
	
	\begin{equation}
	\max_{\boldsymbol{\lambda(x)}} \sum_{\forall x \in X} \lambda(x)
	\end{equation}
	
	s.t. $\forall x \in X$:
	
	\begin{align*}
	\lambda(x_{i})-\alpha\sum_{j=1}^{N}p_{i,j}(u_{1})\lambda(x_{j}) &\leq g(x_{i},u_{1}) \\
	&\vdots\\
	\lambda(x_{i})-\alpha\sum_{j=1}^{N}p_{i,j}(u_{P})\lambda(x_{j}) &\leq g(x_{i},u_{P})
	\end{align*}
	
	Hence, there are $N\cdot P$ constraints in the LP. This form shows that the expected value of the random next state cost is a linear combination of all the costs in $\lambda(S)$, so the constraints are linear.\\
	
	One can rewrite this in the form where the random perturbation $w$ is a state and with the state transition function $f$. Let $W$ be the set of perturbations, indexed by $1<j<M$ . The probability $p_{i,j}(u_{p})$ is defined as:
	
	\begin{align*} 
	p_{i,j}(u_{p})&= P(f(x_{i},u_{p},w_{j}),w_{k} | x_{i},u_{p},w_{j}) \\ 
	&= P(w_{k} | f(x_{i},u_{p},w_{j}),x_{i},u_{p},w_{j})\cdot P(f(x_{i},u_{p},w_{j})| x_{i},u_{p},w_{j})
	\end{align*}
	
	Since the following state for a given combination $x_{i},u_{p},w_{j}$ can only either be feasible or not:
	
	\begin{displaymath}
	P(f(x_{i},u_{p},w_{j})| x_{i},u_{p},w_{j}) = 
	\left\{
	\begin{array}{l}
	1\\
	0
	\end{array}
	\right.
	\end{displaymath}
	
	One may assume it is feasible; otherwise, its term will not be part of the expected value.
	
	\begin{displaymath} 
	p_{i,j}(u_{p})=P(w_{k} | f(x_{i},u_{p},w_{j}),x_{i},u_{p},w_{j})
	\end{displaymath}
	
	%OLD
	%For now, assume the perturbations are IID, so:
	
	%\begin{displaymath} 
	%	p_{i,j}(u_{p})\approx P(w_{k} | f(x_{i},u_{p},w_{j}),x_{i},u_{p})
	%\end{displaymath}\\
	
	%As well, the next perturbation is assumed to only be dependent on the current %state through the next state, so:
	
	%\begin{displaymath} 
	%	p_{i,j}(u_{p})\approx P(w_{k} | f(x_{i},u_{p},w_{j}))
	%\end{displaymath}\\
	
	This can be written just in terms of the current state as:
	\begin{displaymath} 
	p_{i,j}(u_{p})=P(w_{k} | x_{i},u_{p},w_{j})
	\end{displaymath}
	
	%With these assumptions, and ...
	With $w$ being a part of the state, the above LP can be expressed as:
	
	\begin{equation} \label{eq:prelimLP}
	\min_{\boldsymbol{\lambda(x,w)}} \sum_{\forall x \in X,\forall w \in W} -\lambda(x,w)
	\end{equation}
	
	s.t. $\forall x \in X,\forall w \in W$:
	
	%OLD
	%\begin{align*}
	%\lambda(x_{i},w_{j})-\alpha\sum_{k=1}^{M}P(w_{k} | %f(x_{i},u_{1},w_{j}))\lambda(f(x_{i},u_{1},w_{j}),w_{k}) &\leq %g(x_{i},u_{1},w_{j}) \\
	%&\vdots\\
	%\lambda(x_{i},w_{j})-\alpha\sum_{k=1}^{M}P(w_{k} | %f(x_{i},u_{P},w_{j}))\lambda(f(x_{i},u_{P},w_{j}),w_{k}) &\leq %g(x_{i},u_{P},w_{j})
	%\end{align*}
	
	\begin{align*}
	\lambda(x_{i},w_{j})-\alpha\sum_{k=1}^{M}P(w_{k} | x_{i},u_{1},w_{j})\lambda(f(x_{i},u_{1},w_{j}),w_{k}) &\leq g(x_{i},u_{1},w_{j}) \\
	&\vdots\\
	\lambda(x_{i},w_{j})-\alpha\sum_{k=1}^{M}P(w_{k} | x_{i},u_{P},w_{j})\lambda(f(x_{i},u_{P},w_{j}),w_{k}) &\leq g(x_{i},u_{P},w_{j})
	\end{align*}
	
	This is the LP to be solved to determine the converged cost in IHDP. In the context of the hybrid storage problem, $x=(E_{1},E_{2})$ and $w=L$.
	
	The problem can be expressed most succinctly in matrix form as:
	
	\begin{equation}
	\min_{\boldsymbol{\lambda_{t}}} -\boldsymbol{1}_{1\times NM} \boldsymbol{\lambda_{t}}
	\end{equation}
	
	s.t.
	
	\begin{displaymath} 
	\boldsymbol{\lambda_{t}}-\alpha P\boldsymbol{\lambda_{t+1}} \preceq \boldsymbol{g}
	\end{displaymath}
	
	where:
	
	\begin{itemize}
		\item \boldsymbol{$\lambda_{t}$} contains P repeated vectors $\lambda(x,w)$. $dim(\lambda_{t})=NMP\times 1$
		
		\item \boldsymbol{$\lambda_{t+1}$} contains costs for all combinations of next states and next demands. $dim(\lambda_{t+1})=NM^{2}P\times 1$
		
		\item matrix $P$ is a block diagonal  matrix containing the conditioned probabilities of next state demand. $dim(P)=NMP\times NM^{2}P$
		
		\item \boldsymbol{$g$} contains stage cost for each combination of state, control, and current demand. $dim(g)=NMP\times 1$
		
	\end{itemize}
	
	This form cannot yet be solved because there are two cost vectors. One may reduce this to a single decision variable by relating the next state cost to the current state cost.
	
	Firstly, one can consider just the subset of costs associated with the $p^{th}$ control. This is done because each of the sets of linear inequalities in \eqref{eq:prelimLP} ($\forall i,j$) is of the same form, so the full matrix formulation will consist just of the concatenation of the individual matrix equations. For the $p^{th}$ control, the $p^{th}$ set of linear inequalities can be expressed in the form:
	\begin{gather}
	\begin{bmatrix}
	& \lambda_{1}\\
	& \vdotswithin\\
	& \\
	& \lambda_{MN}
	\end{bmatrix}
	-
	\alpha
	\begin{bmatrix}
	P(w_{1}|z_{a}) \hdots  P(w_{M}|z_{a})  \\
	& P(w_{1}|z_{b}) \hdots  P(w_{M}|z_{b}) \\
	& & \ddots
	\end{bmatrix}
	\begin{bmatrix}
	& \lambda_{a,1}\\
	& \vdotswithin\\
	& \\
	& \lambda_{a,M}\\
	& \lambda_{b,1}\\
	& \vdotswithin\\
	& \\
	& \lambda_{b,M}\\
	& \vdotswithin\\
	\end{bmatrix}
	\end{gather}
	
	Here, $z=(x,u,w)$ and the subscripts a and b are used to denote the indices of next possible states (i,j). (Linear indexing is used.)
	Critically, $1\leq a \leq MN $ and $1\leq b \leq MN $, meaning that the next states lie in the same state space as the current states. This implies that it is possible to bijectively map the next state costs $\lambda_{t+1}$ to the current state costs $\lambda_{t}$. A sample rearrangement matrix equation could then take the following form:
	\begin{gather}
	\begin{bmatrix}
	& \lambda_{a,1}\\
	& \vdotswithin\\
	& \\
	& \lambda_{a,M}\\
	& \lambda_{b,1}\\
	& \vdotswithin\\
	& \\
	& \lambda_{b,M}\\
	& \vdotswithin\\
	\end{bmatrix}
	=
	\begin{bmatrix}
	\bold{0} & \hdots & \bold{0} & \bold{I} & \bold{0} & \hdots & \bold{0} \\
	\bold{0} & \bold{I} & \bold{0} & \bold{0} & \hdots & \hdots & \bold{0} \\
	& & & \vdotswithin\\
	\end{bmatrix}
	\begin{bmatrix}
	& \lambda_{1}\\
	& \vdotswithin\\
	& \\
	& \lambda_{1}\\
	& \lambda_{b}\\
	& \vdotswithin\\
	& \\
	& \lambda_{b}\\
	& \vdotswithin\\
	& \\
	& \lambda_{a}\\
	& \vdotswithin\\
	& \\
	& \lambda_{a}\\
	& \vdotswithin\\
	\end{bmatrix}
	\end{gather}
	
	One can define $F$ to be the centre (rearrangement) matrix. Furthermore, the vector on the right can be formed from the current cost vector $\lambda_{t}$ as:
	\begin{gather}
	\begin{bmatrix}
	& \lambda_{1}\\
	& \vdotswithin\\
	& \\
	& \lambda_{1}\\
	& \lambda_{b}\\
	& \vdotswithin\\
	& \\
	& \lambda_{b}\\
	& \vdotswithin\\
	& \\
	& \lambda_{a}\\
	& \vdotswithin\\
	& \\
	& \lambda_{a}\\
	& \vdotswithin\\
	\end{bmatrix}
	=
	\begin{bmatrix}
	\bold{1}  \\
	& \bold{1} \\
	& & \ddots
	\end{bmatrix}
	\begin{bmatrix}
	& \lambda_{1}\\
	& \vdotswithin\\
	& \\
	& \lambda_{MN}
	\end{bmatrix}
	\end{gather}
	
	The central matrix can be termed $Id$.
	
	Hence, the system in \label{eq:prelimLP} can be written as the following LP:
	
	\begin{equation}
	\min_{\boldsymbol{\lambda_{t}}} -\boldsymbol{1}_{1\times NM} \boldsymbol{\lambda_{t}}
	\end{equation}
	
	s.t.
	
	\begin{displaymath} 
	(I-\alpha PFId)\boldsymbol{\lambda_{t}} \preceq \boldsymbol{g}
	\end{displaymath}
	
	\section{Dual LP Problem}
	The above LP does not allow for determining the optimal policy upon convergence, so it is necessary to reformulate it as a minimization over some function of the optimal policy. This can be done by taking the dual of the LP, as shown in \cite{4220813}.
	
	The dual of the LP can be expressed as:
	
	\begin{equation}
	\min_{\boldsymbol{d}} -\boldsymbol{g}^{T} \boldsymbol{d}
	\end{equation}
	
	s.t.
	
	\begin{displaymath} 
	(\Xi-\alpha P^{T})\boldsymbol{d} + \boldsymbol{1} = \boldsymbol{0}
	\end{displaymath}
	
	\begin{displaymath} 
	\boldsymbol{d} \succeq \boldsymbol{0}
	\end{displaymath}
	
	where $dim(\boldsymbol{d})=NMP$.\\
	
	Then, define $\pi$ as a matrix containing in its entries the probabilities of each action $u$ conditioned on state $(x,w)$. The matrix for the optimal sequence, $\pi^{*}$, can be determined as follows, according to \cite{4220813}:
	
	\begin{equation}
	\pi^{*}(x,u)=\frac{d^{*}(x,u)}{\sum_{\forall u \in U}d^{*}(x,u)}
	\end{equation}
	
	Note that when the vertices of the constraints are used, the matrix $\pi$ is purely deterministic\cite{MDPs}. (This is to say that it only indicates optimal controls with certainty.)
	
	
	
	
	\section{Post-Decision Variables}
	Let $y_{1}=f_{1}(E_{1},D_{1})$ and $y_{2}=f_{2}(E_{2},D_{1},D_{2},L)$, where $f_{1}(\cdot)$ and $f_{2}(\cdot)$ are the state equations, respectively:
	
	\begin{displaymath}
	E_{1}(t+1)=\beta_{1}E_{1}(t)+\left[-\frac{1}{\alpha_{1}^{D}}D_{1}(t)\right]
	\end{displaymath}
	
	\begin{displaymath}
	E_{2}(t+1)=\beta_{2}E_{2}(t)+\left[\alpha_{2}^{C}[D_{1}(t)+D_{2}(t)-L(t)]-\frac{1}{\alpha_{2}^{D}}D_{2}(t)\right]
	\end{displaymath}
	
	We can use this to re-write the FHDP problem in \eqref{eq:FHDP_eq} as:
	
	\begin{multline}
	J_{t}[E_{1}(t),E_{2}(t),L(t)] =
	\frac{\alpha_{2}^{C}+\alpha_{2}^{D} -2}{\alpha_{2}^{C}-\frac{1}{\alpha_{1}^{D}}}(\beta_{2}E_{2}(t)+\beta_{1}\alpha_{2}^{C}\alpha_{1}^{D}E_{1}(t))\\
	+\min_{y_{1},y_{2},D_{1}}
	(1-\alpha_{1}^{D})D_{1}(t) + K[D_{1}(t)]^{2} +(1-\alpha_{2}^{C})[D_{1}(t)-L(t)]\\
	-\frac{\alpha_{2}^{C}+\alpha_{2}^{D} -2}{\alpha_{2}^{C}-\frac{1}{\alpha_{1}^{D}}}\left[y_{2}(t)+\alpha_{2}^{C}\alpha_{1}^{D}y_{1}(t)+\alpha_{2}^{C}L(t)\right]\\
	+\mathop{\E}_{L(t+1)}\{J_{t+1}[y_{1}(t),y_{2}(t),L(t+1)]\}
	\end{multline}
	s.t. $\forall j \in \{1,2\}$
	\begin{itemize}
		\item Bounds on stored energy: 
		\begin{math}E_{j}^{min}\leq E_{j}(t)\leq E_{j}^{max}\end{math}
		\item Bounds on next state energy (post-decision variable):
		\begin{math}E_{j}^{min}\leq y_{j}(t)\leq E_{j}^{max}\end{math}
		\item Bounds on discharging:
		\begin{math}0\leq D_{j}(t)\leq D_{j}^{max}\end{math}
		\item Battery state equation:
		\begin{math}y_{1}=f_{1}(E_{1},D_{1})\end{math}
		\item Supercapacitor state equation:
		\begin{math}y_{2}=f_{2}(E_{2},D_{1},D_{2},L)\end{math}
	\end{itemize}
	
	The minimization is of the form:
	\begin{displaymath}
	\min_{y_{1},y_{2},D_{1}} f(y_{1}(t),y_{2}(t),D_{1}(t),L(t)) + \mathop{\E}_{L(t+1)}\{J_{t+1}[y_{1}(t),y_{2}(t),L(t+1)]\}
	\end{displaymath}
	
	Since the state ($E_{j}$) is not in the minimization, one would hope to ensure that its bounds are made redundant simply by only evaluating costs for states in the range \begin{math}E_{j}^{min}\leq E_{j}(t)\leq E_{j}^{max}\end{math}.
	
	This could be done by ensuring that the constraints on $y_{j}$ and $D_{j}$ would automatically correspond to states within this range.
	
	Considering the constraints on $y_{j}$, these automatically set constraints on $E_{j}$ and $D_{j}$:
	\begin{displaymath}E_{1}^{min}\leq f_{1}(E_{1},D_{1})\leq E_{1}^{max}\end{displaymath}
	\begin{displaymath}E_{2}^{min}\leq f_{2}(E_{2},D_{1},D_{2},L)\leq E_{2}^{max}\end{displaymath}
	
	Substituting for $f_{1}(\cdot)$ and $f_{2}(\cdot)$:
	\begin{displaymath}E_{1}^{min}\leq \beta_{1}E_{1}(t)+\left[-\frac{1}{\alpha_{1}^{D}}D_{1}(t)\right]\leq E_{1}^{max}\end{displaymath}
	\begin{displaymath}E_{2}^{min}\leq \beta_{2}E_{2}(t)+\left[\alpha_{2}^{C}[D_{1}(t)+D_{2}(t)-L(t)]-\frac{1}{\alpha_{2}^{D}}D_{2}(t)\right]\leq E_{2}^{max}\end{displaymath}
	
	In the first constraint, one sees that it is possible to bound the control $D_{1}$ in by the variable bound of the state $E_{1}$.
	
	However, the control $D_{2}$ in the second constraint is not \textit{independently} defined for this minimization, and so the constraint on state $E_{2}$ cannot be made redundant.
	
	Therefore, Post-Decision Variables \textit{cannot} be used to reduce the number of discretized states for this problem.
	
	
	
	\section{Miscellaneous Notes}
	Single storage:
	\begin{multline}
	J_{t}[E_{2}(t),L(t)] = \min_{D_{1},D_{2}}
	(1-\alpha_{1}^{D})D_{1}(t) 
	+ K[D_{1}(t)]^{2}
	+(1-\alpha_{2}^{C})[D_{1}(t)-L(t)]\\
	+\mathop{\E}_{L(t+1)} \{J_{t+1}[f_{2}(E_{2}(t),D_{1}(t)),L(t+1)]\}
	\end{multline}
	Iteration costs:
	\begin{itemize}
		\item $N^{th}$ stage: \\
		\begin{equation}
		J_{N-1}[x(N-1),w(N-1)]=0
		\end{equation}
		
		\item $(N-1)^{th}$ stage: \\
		\begin{multline}
		J_{N-2}[x(N-2),w(N-2)]=\min_{u} g(x(N-2),u(N-2),w(N-2))+\\ \mathop{\E}_{w(N-1)}\{0\}
		\end{multline}
		
		\item $(N-2)^{th}$ stage: \\
		\begin{multline}
		J_{N-3}[x(N-3),w(N-3)]=\min_{u} g(x(N-3),u(N-3),w(N-3))\\
		+ \mathop{\E}_{w(N-1)} \{
		g( f(x(N-3),u(N-3),w(N-3)) ,u(N-2),w(N-2)) + 0\}
		\end{multline}
		
		\item Initial stage: \\	
		\begin{multline}
		J=J_{0}[x(0),w(0)]=\min_{u} g(x(0),u(0),w(0))\\
		+\left[\sum_{t=1}^{N-1}\mathop{\E}_{w(t)} \{
		g( f(x(t-1),u(t-1),w(t-1)) ,u(t),w(t))
		\}\right]
		\end{multline}
	\end{itemize}
	
	Note that here, $f(\cdot)$ is a general recursive state equation.
	
	
	%\bibliographystyle{plain}
	%\bibliography{references}
	
	%\printbibliography
\end{document}