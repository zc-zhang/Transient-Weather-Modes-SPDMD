Visualization complex valued modes
abs(ϕ) (Magnitude):

Use: It shows the overall amplitude of the Koopman modes and is useful for identifying dominant modes.
Transient Behavior: While it can indicate where the mode is most active, it might not fully capture transient dynamics since it ignores phase information and the direction of oscillations.
Re(ϕ) (Real Part):

Use: Captures the component of the mode that corresponds to steady-state behavior or slowly varying oscillations.
Transient Behavior: The real part can highlight how the mode contributes to persistent or slowly decaying transients.
Im(ϕ) (Imaginary Part):

Use: Represents the oscillatory component, often associated with cycles or rotations in the state space.
Transient Behavior: If the system exhibits significant oscillatory transients, the imaginary part can be quite informative.
∠(ϕ) (Phase):

Use: Reflects the phase relationship and can show how different modes synchronize or interact.
Transient Behavior: Phase information is crucial in understanding the timing of oscillations, especially in coupled or interacting systems.

% load 'scalehandle240805b.mat'
% run 'vorflow.m' to get magnitude of vorticity (scale filed) data, i.e., 'Vfull.mat' with size 3880*121

% load 'Vfull.mat'
% run 'spdmd_vorscale.m' to perform sparsity-promoting dynamic mode decomposition (spDMD) proposed by MR Jovanovic, Physics of Fluids'2014.

% Plot figures:
% run 'plot_spdmd_vor.m'
% run 'plot_STModes_spdmd.m' 

% when plotting it need the color-shadow
% use 'shadedErrorBar.m' to show the shadow
 
