# medical trial Monte Carlo simulation
The primary purpose of this project is documentation; only a minimum level of comfort is provided. Is is used to estimate the probability of success of a specific clinical trial (namely, the Sellas Live Sciences REAGL phase 3 trial) using Monte Carlo simulations. The simulation can be performed using the program SLSMonteCarlo.cpp. The averaged point of success (PoS) is displayed on the command line, and partial results are written to the file simulation_results.csv. These partial results can be visualized using SlsMonteCarloEvaluation.m. For more details on the background see my blog post on [Monte Carlo simulations for the REGAL clinical trial](https://vola-corvidae.com/artikel/monte_carlo_simulations_sls.html).

# Installation and usage
Integrate SLSMonteCarlo.cpp into your preferred build system. The file depends on the Boost library. Adapt the parameters as described below. Compile and run. Depending on your machine and choice of parameters runtime can be around one hour. Read the main result from command line. Use the octave script if visualize the marginalized mOS and shape grids.

# Parameters
Before you compile review the parameters at line 34 and following.<br>
*constexpr bool calcPower<br>
If set filters are not applied this is equivalent to calculate the initial power of the trial.<br>
<br>
constexpr bool shapeGrid;<br>
When this option is enabled, the simulation performs a mesh calculation using a certain range of shape parameters. This increases the runtime by a factor of 49.<br>
<br>
constexpr double shapeDefault = 0.5;<br>
Only used if shapeGrid is set to false. Sets the value of the shape parameter that is used for all simulations.<br>
<br>
const size_t numSimulations = calcPower ? 1000 : 10000;<br>
How many simulations are performed for each parameter combination.<br>
<br>
constexpr double maxTime = 60.0; // max 5 years (60 months)<br>
Maximum survival time taken into account.<br>
<br>
const size_t maxPatientsPerArm = 63;<br>
self-explanatory<br>
<br>
const double minCondPowerFutility = 0.1;<br>
What is the minimum required conditional power in the interim analysis to pass the futility filters<br>
<br>
constexpr double timeAtIa1 = 8; // after end of enrollment<br>
Date of the first interim analysis in months after the end of the recruitment phase.<br>
<br>
constexpr double timeAtIa2 = 16.0; // after end of enrollment<br>
Date of the second interim analysis in months after the end of the recruitment phase.<br>
<br>
constexpr double timeCurrent = 20.0; // after end of enrollment<br>
Time point up to which the 80th event has not occurred.<br>
<br>
constexpr double timeEnrollment = 14; // months<br>
Duration of the recruitment phase in months.<br>
<br>
constexpr double graceFilter1 = 3.0; // months<br>
How much the date of the first interim analysis may deviate from the actual date in months.<br>
<br>
constexpr double medianFollowUpTimeIa1 = 13.5; // months<br>
The median follow up time at first interim analysis in months. 13.5 month are taken from Sellas Live Sciences [press release](https://ir.sellaslifesciences.com/news/News-Details/2025/SELLAS-Life-Sciences-Announces-Positive-Outcome-of-Interim-Analysis-for-its-Pivotal-Phase-3-REGAL-Trial-of-GPS-in-Acute-Myeloid-Leukemia/default.aspx)<br>
constexpr double graceMedianFollowUpTime = 3.0;<br>
How much the median follow up time at the first interim analysis may deviate from the actual date in months.
<br>
const size_t eventsIaAll = 60;<br>
Number of events that trigger the first interim analysis.<br>
<br>
const size_t eventsFaAll = 80;<br>
Number of events that trigger thefinal analysis.<br>
<br>
const double pLimitFa = 0.047;<br>
Expected alpha spend on final analysis.<br>
<br>
const double pLimitIa = 0.003;<br>
Expected alpha spend on first interim analysis.<br>
<br>
const double median_control_prior = 10.0; // mean overall survival (mOS) (months)<br>
Center for the prior of the mean overall survival  of the control (BAT) arm in month.<br>
<br>
const double std_control_prior = 2.5;   // standard deviation (control)<br>
Standard deviation for the prior of the mean overall survival  of the control (BAT) arm in month.<br>
<br>
const double median_treatment_prior = 18.0; // mOS<br>
Center for the prior of the mean overall survival  of the treatment (GPS) arm in month.<br>
<br>
const double std_treatment_prior = 6.0;   // standard deviation (treatment)<br>
Standard deviation for the prior of the mean overall survival  of the treatment (GPS) arm in month.<br>
