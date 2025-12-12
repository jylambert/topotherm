Methodology
=============



Linearization of heat losses and pipe investment costs
--------------------------------------------------------
 The maximal mass flow :math:`\dot{m}`  with the corresponding velocities :math:`v` of each considered piping diameter :math:`d` must be determined to linearize the maximal thermal power flow in a district heating pipe with length l. Subsequently, the heat losses of the district heating network pipes can be determined for each diameter. The relationship between the investment costs and heat losses of a pipe with regard to the maximal thermal power flow is modeled with a linear regression.

Initial Assumptions
^^^^^^^^^^^^^^^^^^^^

 * The considered heat carrier medium is liquid water and, therefore, incompressible
 * Temperature changes are assumed to be smaller than :math:`\Delta`T = 40 °C, neglecting temperature dependencies of the fluid properties (e.g.,density :math:`\rho`, heat capacity and dynamic viscosity :math:`\mu`).

Calculations
^^^^^^^^^^^^^
The following equations are used to form a full non-linear thermo-hydraulic model of a district heating pipe.

Calculating the pressure loss between the nodes :math:`i` and :math:`j`:

.. math::
  
  \Delta p_{ij} = f_{ij} \cdot \frac{\text{l}_{ij}}{d_{ij}} \cdot \frac{v_{ij}^2}{2} \cdot \rho
   
Calculating the Reynolds Number:

.. math:: 
  
  Re_{ij} = \frac{\rho \cdot v_{ij} \cdot d_{ij}}{\mu} 

Calculate the friction factor for a turbulent flow (Re > 2300) with pipe roughness :math:`\epsilon`

.. math:: 

  f_{ij} = \left[ -1.8 \cdot \log \left( \left(\frac{\varepsilon}{3.7 \cdot d_{ij}}\right)^{1.11} + \frac{6.9}{Re_{ij}} \right) \right]^{-2}
    
An iterative calculation (using SciPy) is performed to determine the maximum velocity in a pipe with a given inner diameter. The specific pressure drop per meter pipe length should range from 70 Pa/m to 350 Pa/m and the initial velocity is 0.5 m/s. 

.. math:: 
  
  v_{max} = \sqrt{\frac{2 \cdot \Delta p_{max} \cdot d}{f \cdot \rho}}
  
After calculating the maximum speed the maximum flow rate can be calculated using :math:`\dot{m} = \rho \cdot v \cdot A` with the pipe's cross-sectional area A.

To model the thermal behavior of an insulated pipe buried underground the temperature difference :math:`\Theta` between the water temperature in the pipe and the outside temperature is introduced:

.. math:: 
 
 \Theta_{ij} = \Theta_{i} \cdot \exp\left(\frac{{-\text{l}_{ij}}{-l_{ij}}}{c_p \cdot \dot{m}_{ij} \cdot R_{ij}}\right)


The combined thermal resistance :math:`R_{ij}` of the pipe and soil per unit length is calculated with the ratio :math:`r` between outer and inner diameter. The buried depth of a pipe :math:`h`, as well as the combined thermal resistance.

.. math:: 
  
 R_{ij} = \frac{\ln{\frac{4h}{r \cdot d_{ij}}}}{2\pi \cdot k_{g}} + \frac{\ln{r}}{2\pi \cdot k_{insul}}


Linear regression
^^^^^^^^^^^^^^^^^^

To represent the piping diameters in a MILP, they are simplified to continuous and linear. The linearization has the form of :math:`y = a \cdot x + b`, where :math:`y` is the dependent variable, :math:`x` is the independent variable, and the factors :math:`a` and :math:`b` are the regression coefficients.

The thermal losses are calculated through the temperature difference between the flowing water and the pipe calculated earlier using the maximal mass flow rate, and the total investment costs can be assigned by the user.

.. figure:: ../img/regression.png

Mixed-integer formulations for district heating network design
---------------------------------------------------------------

Graphical representation
^^^^^^^^^^^^^^^^^^^^^^^^^

A graph representation of a district heating network is used to represent pipes, internal nodes, supply, and sinks. The pipes correspond to the graph's edges and the network's junctions to the nodes. The network consists of a feed and return network with edges of opposite directions. This superstructure contains all possible connections and pathways from the heat source to the sinks.

.. figure:: ../img/graph.png

Single Time Step Formulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The bidirectionality of the network (two edges per pipe) is explicitly considered and modeled in the single time step formulation of topotherm. In order to allow flows in the opposite direction, every potential pipe :math:`ij` is also modeled in the direction :math:`ji`.

* The heat balance of the pipe

.. math:: 
  
  \dot{Q}_{ij,0} - \dot{Q}_{ij,1} -\left(a_{therm} \cdot \dot{Q}_{ij,0} + b_{therm} \cdot \lambda_{ij} \right) \cdot \text{l}_{ij}  = 0  \qquad \forall ij \in \mathcal{A}_{int}

Note: thermal losses are determined by the linear regression coefficients :math:`a_{therm}` and :math:`a_{therm}` while the binary variable :math:`\lambda_{ij}` represents the flow direction of the considered pipe.

* A big-M constraint is formulated to enforce zero thermal power flow if the direction :math:`ij` is not used.

.. math:: 

 \dot{Q}_{ij,0} \leq \dot{Q}_{max,cons} \cdot \lambda_{ij} \qquad \forall ij \in {\mathcal{A}_{int}}


* Each consumer connection to the district heating grid is modeled as unidirectional, and thus no heat feed-in from a consumer is possible. Consequently, the node energy conservation is formulated as follows:

.. math:: 

  \dot{Q}_{ni,0} - \dot{Q}_{jn,1} - \dot{Q}_{c,n} + \dot{Q}_{p,n} = 0 \qquad \forall n \in \mathcal{N}, \; (ni, jn) \in  {\mathcal{A}_{int}} 

* Enforcement of unidirectional use of a pipe:

.. math:: 

  \lambda_{ij} + \lambda_{ji} \leq 1 \qquad \forall ij \in {\mathcal{A}_{int}} \; 

The annuity method distributes investment costs of pipes or heat sources over the defined life span :math:`n_{years}` with an interest rate :math:`w`.

.. math:: 
 
  an = \frac{(1+w)^{n_{years}} \cdot w}{(1+w)^{n_{years}} -1}
 
Finally, the objective function minimizes the district heating network's total investment and operational costs. The investment costs are determined by the linear regression factors :math:`a_{cost}` and :math:`a_{cost}`. By introducing full load hours flh, the investment and operational costs are weighted.

.. math:: 

 \text{min} \; \bigg\{& \sum_{p \in {\mathcal{N}_p}} {\dot{Q}_{\text{inst},p}} \cdot {c_{\text{inv},p}} \cdot an_p  + \sum_{p \in {\mathcal{N}_p}} \dot{Q}_{p} \cdot {c_{\text{fuel},p}} \cdot {\text{flh}}  \\ & + \sum_{ij \in \mathcal{A}_{int}} \left(a_{cost} \cdot (\dot{Q}_{ij,0}+\dot{Q}_{ji,0}) + b_{cost} \cdot (\lambda_{ij}+\lambda_{ji}) \right) \cdot {\text{l}_{ij}} \cdot an_{pipe} \bigg\}

Multiple Time Step Formulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to optimize district heating systems with multiple time steps and changing flow directions, the sts formulation needs to be expanded. Not only the bidirectional operation of the pipe has to be modeled, but also a binary variable has to be incorporated to model the decision if a pipe is built independent of the flow direction. As with the sts formulation, only equations in direction :math:`ij` are provided below and :math:`t` refers to the set of all considered hourly timesteps.

The constraints of mts are as follows:

* The heating power balance of each pipe

.. math:: 

  \dot{Q}_{ij,0}[t] - \dot{Q}_{ij,1}[t] - {\dot{Q}_{\text{loss},ij}[t]} = 0  \qquad \forall ij \in {\mathcal{A}_{int}}, \forall t \in \mathcal{T} 

* The thermal losses for pipe :math:`ij` need to be modeled as an independent variable and cannot be incorporated into the equation directly.

.. math:: 

  &{\dot{Q}_{\text{loss},ij}[t]} - \dot{Q}_{max,cons} \cdot \lambda_{ij}[t] \leq 0 & \forall ij \in {\mathcal{A}_{int}}, \forall t \in \mathcal{T} \label{eq:loss_1}\\ &\left(a_{therm} \cdot {\dot{Q}_{\text{cap},ij}} + b_{therm} \cdot \lambda_{ij}[t]\right) \cdot {\text{l}_{ij}} - {\dot{Q}_{\text{loss},ij}[t]} - \dot{Q}_{max,cons} \cdot \left(1-\lambda_{ij}[t]\right) \leq 0 & \forall ij \in {\mathcal{A}_{int}}, \forall t \in \mathcal{T}

If the direction 𝑖𝑗 is not used,the first equation ensures that the heat losses equal 0. The second equation enforces the heat loss calculation with the maximal built thermal capacity ̇ 𝑄max according to the flow direction in the pipe.

* Zero thermal flow if the direction :math:`ij` of a pipe is not used.

.. math:: 

  \dot{Q}_{ij,0}[t] \leq \dot{Q}_{max,cons} \cdot \lambda_{ij}[t] \qquad \forall ij \in {\mathcal{A}_{int}}, \forall t \in \mathcal{T}

* The maximal thermal power inflow ̇at each time step is limited to the maximal thermal capacity of a pipe and to the binary decision.

.. math:: 

  &\dot{Q}_{ij,0}[t] \leq {\dot{Q}_{\text{cap},ij}} \qquad & \forall ij \in {\mathcal{A}_{int}}, \forall t \in \mathcal{T} \\ & {\dot{Q}_{\text{cap},ij}} \leq \dot{Q}_{max,cons} \cdot {\lambda_{\text{built},ij}} \qquad & \forall ij \in {\mathcal{A}_{int}}

Objective function:

.. math::

  \text{min} \;  \bigg\{& \sum_{p \in \mathcal{N}_{p}} \dot{Q}_{\text{inst},p} \cdot c_{\text{inv},p} \cdot an_p  + \sum_{t \in \mathcal{T}} \sum_{p \in \mathcal{N}_{p}} \dot{Q}_{p}[t] \cdot c_{\text{fuel},p} \cdot \text{flh}   \\ & + \sum_{ij \in \mathcal{A}_{int}} \left(a_{cost} \cdot \dot{Q}_{\text{cap},ij} + b_{cost} \cdot \lambda_{\text{built},p} \right) \cdot \text{l}_{ij} \cdot an_{pipe} \bigg\}


The forced expansion is modified for an economic expansion by eliminating the constraints of all sink connections and modifying the objective function with economic indicators, such as revenue maximization.
