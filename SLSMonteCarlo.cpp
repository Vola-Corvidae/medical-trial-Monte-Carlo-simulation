#include <boost/random/mersenne_twister.hpp>
#include <boost/random/weibull_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>


#include <stdio.h>
#include <cmath>
#include <numbers>
#include <numeric>
#include <algorithm>
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>
#include <variant>
#include <fstream>

enum Arm
{
   Bat,
   Gps
};
constexpr size_t arms = 2;

struct AtCalendarTime { double time; };
struct AtEventCount { size_t count; };
using AnalysisTrigger = std::variant<AtCalendarTime, AtEventCount>;

// parapeters of simulation
constexpr bool calcPower = false;
constexpr bool shapeGrid = true;
constexpr double shapeDefault = 0.5;
const size_t numSimulations = calcPower ? 1000 : 10000;
constexpr double interval = 0.25; // 0.25 months
constexpr double maxTime = 60.0; // max 5 years (60 months)
constexpr size_t numIntervals = ((size_t)(maxTime / interval));
const size_t maxPatientsPerArm = 63;
const double minCondPowerFutility = 0.1;
constexpr double timeAtIa1 = 8; // after end of enrollment
constexpr double timeAtIa2 = 16.0; // after end of enrollment
constexpr double timeCurrent = 20.0; // after end of enrollment
constexpr double timeEnrollment = 14; // months
constexpr double graceFilter1 = 3.0; // months
constexpr double medianFollowUpTimeIa1 = 13.5; // months
constexpr double graceMedianFollowUpTime = 3.0;

const size_t eventsIaAll = 60;
const size_t eventsFaAll = 80;
const double pLimitFa = 0.047;
const double pLimitIa = 0.003;
 
const double median_control_prior = 10.0; // mean overall survival (mOS) (months)
const double std_control_prior = 2.5;   // standard deviation (control)
const double median_treatment_prior = 18.0; // mOS
const double std_treatment_prior = 6.0;   // standard deviation (treatment)

boost::random::mt19937 rng(static_cast<unsigned int>(time(NULL)));

double random_uniform(double lower, double upper) {
   static boost::random::uniform_real_distribution<double> dist(lower, upper);
   return dist(rng);
}

double weibull_random(double scale, double shape) {
   boost::random::weibull_distribution<double> dist(shape, scale);
   return dist(rng);
}

double lognormal_pdf(double x, double mu, double sigma) {
   if (x <= 0) return 0.0;
   double log_x = log(x);
   double exponent = -pow(log_x - mu, 2) / (2.0 * sigma * sigma);
   return (1.0 / (x * sigma * sqrt(2.0 * std::numbers::pi))) * exp(exponent);
}

// Convert expectation value and standard deviation to mu and sigma.
void lognormal_params(double median, double std, double* mu, double* sigma) {
   if (median <= 0 || std <= 0) {
      *mu = 0.0;
      *sigma = 1.0;
      return;
   }
   *sigma = sqrt(log(1.0 + (std * std) / (median * median)));
   *mu = log(median);
}

double chiSquaredToPvalue(double chiSquaredValue, double degreesOfFreedom)
{
   if (degreesOfFreedom <= 0) {
      throw std::invalid_argument("Degrees of freedom must be positive");
   }

   boost::math::chi_squared dist(degreesOfFreedom);
   double cdf_value = boost::math::cdf(dist, chiSquaredValue);
   return 1.0 - cdf_value;
}

std::tuple<std::vector<double>, double, size_t> countEvents(
   const std::array<std::array<double, maxPatientsPerArm>, arms>& survivalTimes,
   const std::array<std::array<double, maxPatientsPerArm>, arms>& enrollmentTimes,
   AnalysisTrigger trigger)
{
   std::vector<double> allEvents;

   for (size_t a = 0; a < arms; ++a) {
      for (size_t i = 0; i < maxPatientsPerArm; ++i) {
         if (!std::isfinite(survivalTimes[a][i])) [[unlikely]] continue;

         double eventCalendarTime = enrollmentTimes[a][i] + survivalTimes[a][i];
         allEvents.push_back(eventCalendarTime);
      }
   }

   if (allEvents.empty()) return { {}, 0.0, 0 };

   std::sort(allEvents.begin(), allEvents.end(),
      [](const double& a, const double& b) { return a < b; });
   double timeOfAnalysis = std::numeric_limits<double>::infinity();
   
   if (std::holds_alternative<AtEventCount>(trigger)) {
      size_t target = std::get<AtEventCount>(trigger).count;
      if (target == 0)
         return { {}, 0.0, 0 };
      if (target > allEvents.size()) {
         timeOfAnalysis = allEvents.back();
      }
      else {
         timeOfAnalysis = allEvents[target - 1];
      }
   }
   else if (std::holds_alternative<AtCalendarTime>(trigger)) {
      timeOfAnalysis = std::get<AtCalendarTime>(trigger).time;
   }
   
   std::vector<double> riskSetTimes;
   riskSetTimes.reserve(allEvents.size());
   for (const auto& ev : allEvents) {
      if (ev > timeOfAnalysis + 1e-10) break;
      riskSetTimes.push_back(ev);
   }
   size_t events = riskSetTimes.size();

   return { riskSetTimes, timeOfAnalysis, events };
}

struct LogRankComponents {
   double O1 = 0.0;
   double E1 = 0.0;
   double V = 0.0;
   size_t totalEvents = 0;
   double timeOfAnalysis;
};

std::array<std::vector<double>, arms> sortedExitTimes(
   const std::array<double, maxPatientsPerArm* arms>& survivalTimes,
   const std::array<double, maxPatientsPerArm* arms>& enrollmentTimes)
{
   constexpr size_t totalPatients = maxPatientsPerArm * arms;
   std::array<std::vector<double>, arms> exitTimes;
   for (int a = 0; a < arms; ++a) exitTimes[a].reserve(maxPatientsPerArm);

   for (size_t i = 0; i < totalPatients; ++i) {
      if (!std::isfinite(survivalTimes[i])) [[unlikely]] continue;
      int arm = (i < maxPatientsPerArm) ? 0 : 1;
      exitTimes[arm].push_back(enrollmentTimes[i] + survivalTimes[i]);
   }

   for (int a = 0; a < arms; ++a)
      std::sort(exitTimes[a].begin(), exitTimes[a].end());
   return exitTimes;
}

LogRankComponents computeLogRankComponents(
   const std::array<double, maxPatientsPerArm* arms>& survivalTimes,
   const std::array<double, maxPatientsPerArm* arms>& enrollmentTimes,
   std::array<std::vector<double>, arms> exitTimes,
   AnalysisTrigger trigger)
{
   LogRankComponents out{};

   double timeOfAnalysis = 0.0;
   std::vector<double> allEvents;
   allEvents.reserve(exitTimes[0].size() + exitTimes[1].size());
   allEvents.insert(allEvents.end(), exitTimes[0].begin(), exitTimes[0].end());
   allEvents.insert(allEvents.end(), exitTimes[1].begin(), exitTimes[1].end());
   std::sort(allEvents.begin(), allEvents.end()); // nur für AtEventCount nötig

   if (std::holds_alternative<AtEventCount>(trigger)) {
      size_t target = std::get<AtEventCount>(trigger).count;
      if (target == 0 || target > allEvents.size())
         timeOfAnalysis = allEvents.back();
      else
         timeOfAnalysis = allEvents[target - 1];
   }
   else {
      timeOfAnalysis = std::get<AtCalendarTime>(trigger).time;
   }

   out.timeOfAnalysis = timeOfAnalysis;

   size_t idx[2] = { 0,0 };
   size_t totalEvents = 0;

   while (idx[0] < exitTimes[0].size() || idx[1] < exitTimes[1].size()) {
      double nextEvent;

      if (idx[0] >= exitTimes[0].size()) nextEvent = exitTimes[1][idx[1]];
      else if (idx[1] >= exitTimes[1].size()) nextEvent = exitTimes[0][idx[0]];
      else nextEvent = std::min(exitTimes[0][idx[0]], exitTimes[1][idx[1]]);

      if (nextEvent > timeOfAnalysis + 1e-10)
         break;


      size_t d[2] = { 0,0 };
      while (idx[0] < exitTimes[0].size() && std::abs(exitTimes[0][idx[0]] - nextEvent) < 1e-10) { ++d[0]; ++idx[0]; }
      while (idx[1] < exitTimes[1].size() && std::abs(exitTimes[1][idx[1]] - nextEvent) < 1e-10) { ++d[1]; ++idx[1]; }

      size_t n[2] = { exitTimes[0].size() - idx[0] + d[0],
                      exitTimes[1].size() - idx[1] + d[1] };

      size_t N = n[0] + n[1];
      size_t D = d[0] + d[1];

      if (N <= 1 || D == 0) continue;

      double e1 = static_cast<double>(n[1]) * D / N; // Arm 1 = Gps
      out.O1 += d[1];
      out.E1 += e1;
      out.V += (static_cast<double>(n[0]) * n[1] * D * (N - D)) / (N * N * (N - 1.0));

      totalEvents += D;
   }

   out.totalEvents = totalEvents;
   return out;
}

double logRankStatistic(const LogRankComponents& lrStat)
{
   if (lrStat.V < 1e-12) 
      return 0.0;
   return (lrStat.O1 - lrStat.E1) * (lrStat.O1 - lrStat.E1) / lrStat.V;
}

double calculateHazardRatio(
   std::array<double, maxPatientsPerArm * arms>& survivalTimes,
   std::array<double, maxPatientsPerArm * arms>& enrollmentTimes,
   std::array<std::vector<double>, arms> exitTimes,
   AnalysisTrigger trigger)
{
   auto C = computeLogRankComponents(survivalTimes, enrollmentTimes, exitTimes, trigger);

   if (C.V < 1e-12 || C.totalEvents == 0)
      return 1.0;

   double logHR = (C.O1 - C.E1) / C.V;
   return std::exp(logHR);
}

double calculateMedianFollowUpTime(
   const std::array<double, maxPatientsPerArm* arms>& survivalTimes,
   const std::array<double, maxPatientsPerArm* arms>& enrollmentTimes,
   AnalysisTrigger trigger)
{
   constexpr size_t totalPatients = maxPatientsPerArm * arms;

   std::array<double, totalPatients> eventCalendarTimes{};
   size_t eventCount = 0;

   for (size_t i = 0; i < totalPatients; ++i) {
      if (!std::isfinite(survivalTimes[i])) continue;
      double eventTime = enrollmentTimes[i] + survivalTimes[i];
      eventCalendarTimes[eventCount++] = eventTime;
   }

   if (eventCount == 0) return 0.0;

   std::sort(eventCalendarTimes.begin(), eventCalendarTimes.begin() + eventCount);

   double timeOfAnalysis = 0.0;
   if (std::holds_alternative<AtEventCount>(trigger)) {
      size_t target = std::get<AtEventCount>(trigger).count;
      if (target == 0 || target > eventCount) {
         timeOfAnalysis = eventCalendarTimes[eventCount - 1];
      }
      else {
         timeOfAnalysis = eventCalendarTimes[target - 1];
      }
   }
   else {
      timeOfAnalysis = std::get<AtCalendarTime>(trigger).time;
   }

   alignas(64) static thread_local std::array<double, totalPatients> followUpBuffer{};
   size_t validCount = 0;

   for (size_t i = 0; i < totalPatients; ++i) {
      if (!std::isfinite(survivalTimes[i])) continue;

      double entry = enrollmentTimes[i];

      // Patient must have been recruited before or at the time of analysis.
      if (entry >= timeOfAnalysis + 1e-10) continue;

      double observedTime = survivalTimes[i];
      double censoredAt = timeOfAnalysis - entry;

      double followUp = std::min(observedTime, censoredAt);
      if (followUp > 0.0) {
         followUpBuffer[validCount++] = followUp;
      }
   }

   if (validCount == 0)
      return 0.0;

   // Calculate the median (partial sorting is sufficient)
   // We only need the median → nth_element is faster than full sort
   size_t mid = validCount / 2;

   std::nth_element(followUpBuffer.begin(), followUpBuffer.begin() + mid, followUpBuffer.begin() + validCount);

   if (validCount % 2 == 1) {
      return followUpBuffer[mid];
   }
   else {
      // Find the second smallest one in the upper half
      double val1 = followUpBuffer[mid - 1];
      double val2 = followUpBuffer[mid];

      // If already sorted, than val1 <= val2
      // Otherwise, find the maximum of the lower half and the minimum of the upper half.
      if (val1 > val2) std::swap(val1, val2);

      // Better: second value set to mid-1 by nth_element
      std::nth_element(followUpBuffer.begin(), followUpBuffer.begin() + mid - 1, followUpBuffer.begin() + validCount);
      val1 = followUpBuffer[mid - 1];

      return (val1 + val2) * 0.5;
   }
}

double calcConditionalPower(const double chiSquare, const double nIa, const double nCont, const double alpha)
{
   const double nF = nCont + nIa;
   const double zAlpha = boost::math::quantile(boost::math::normal_distribution<double>(0.0, 1.0), 1.0 - alpha / 2.0);
   if (chiSquare <= 0 || !std::isfinite(chiSquare)) {
      return 0.0;
   }
   const double zIa = sqrt(chiSquare);
   const double z = (zIa*sqrt(nF / (nF - nCont)) - zAlpha) / sqrt(nF / nCont);
   static boost::math::normal_distribution<double> standardNormalDist(0.0, 1.0);
   return boost::math::cdf(standardNormalDist, z);
}

std::tuple<bool, double> simulateInstance(
   const std::array<double, arms>& medianSurvival,
   const std::array<double, arms>& shape,
   const std::array<size_t, arms>& numPatients
   )
{
   // Initialisiere Patienten
   std::array<double, maxPatientsPerArm * arms> survivalTimes;
   std::array<double, maxPatientsPerArm * arms> enrollmentTimes;
   for (size_t a = 0; a < arms; ++a) {
      double scale = medianSurvival[a] / std::pow(std::log(2), 1.0 / shape[a]);
      for (size_t i = 0; i < numPatients[a]; ++i) {
         size_t idx = a == 0 ? i : (maxPatientsPerArm + i);
         survivalTimes[idx] = weibull_random(scale, shape[a]);
         enrollmentTimes[idx] = random_uniform(0.0, timeEnrollment);
      }
   }

   auto exitTimes = sortedExitTimes(survivalTimes, enrollmentTimes);

   if(!calcPower){
       //   
       // size_t sumEvents = countEvents(intervallsAtIa1);
       // if(sumEvents > eventsIaAll + 10 || sumEvents < eventsIaAll - 10)
       //    return {false, -1.0};
        
        // filter 1 not halted for futility at ia1
        auto lrStatsIa = computeLogRankComponents(survivalTimes, enrollmentTimes, exitTimes, AtEventCount{ eventsIaAll });
        auto chiSquaredIa = logRankStatistic(lrStatsIa);
        if(lrStatsIa.totalEvents != eventsIaAll){
            std::cout << "Error, events Ia: " << lrStatsIa.totalEvents << std::endl;
            abort();
        }
        // chi^2 is symetric, so we have to check for hr > 1 seperatly
        double hrIa = calculateHazardRatio(survivalTimes, enrollmentTimes, exitTimes, AtEventCount{ eventsIaAll });
        if (hrIa > 1.0)
           return { false, -1.0 };
        double condPowerIa = calcConditionalPower(chiSquaredIa, eventsIaAll, eventsFaAll - eventsIaAll, pLimitFa);
        if (condPowerIa < minCondPowerFutility)
           return { false, -1.0 };
        // filter 2 first ia at a specific time after enrollment ended
        if (lrStatsIa.timeOfAnalysis < timeAtIa1 + timeEnrollment - graceFilter1)
           return { false, -1.0 };
        if (lrStatsIa.timeOfAnalysis > timeAtIa1 + timeEnrollment + graceFilter1)
           return { false, -1.0 };
        // filter 3 exactly 60 events at 13.5 months median follow up, (including grace period to keep enough trials)
        auto medianFollowUpTime = calculateMedianFollowUpTime(survivalTimes, enrollmentTimes, AtEventCount{ eventsIaAll });
        if (medianFollowUpTime < medianFollowUpTimeIa1 - graceMedianFollowUpTime )
           return { false, -1.0 };
        if (medianFollowUpTime > medianFollowUpTimeIa1 + graceMedianFollowUpTime)
           return { false, -1.0 };
      
        // filter 4 not halted for efficency at ia1   
        double pIa = chiSquaredToPvalue(chiSquaredIa, 1) / 2.0;
        if(pIa < pLimitIa || pIa > 1.0)
           return { false, -1.0 };
         
        // filter 5 Fa not reached, has to be applied bevor filter 5 otherwise event > 80 will result to crashes
        auto lrStatsCurrent = computeLogRankComponents(survivalTimes, enrollmentTimes, exitTimes, AtCalendarTime{ timeCurrent + timeEnrollment });
        if (lrStatsCurrent.totalEvents >= eventsFaAll)
           return { false, -1.0 };
         
        // filter 6 not halted for futility at ia2
        auto lrStatsIa2 = computeLogRankComponents(survivalTimes, enrollmentTimes, exitTimes, AtCalendarTime{ timeAtIa2 + timeEnrollment });
        auto chiSquaredIa2 = logRankStatistic(lrStatsIa2);
        double hrIa2 = calculateHazardRatio(survivalTimes, enrollmentTimes, exitTimes, AtCalendarTime{ timeAtIa2 + timeEnrollment });
        if (hrIa2 > 1.0)
           return { false, -1.0 };
        condPowerIa = calcConditionalPower(chiSquaredIa2, static_cast<double>(lrStatsIa2.totalEvents), static_cast<double>(eventsFaAll - lrStatsIa2.totalEvents), pLimitFa);
        if (condPowerIa < minCondPowerFutility){
           return { false, -1.0 };
        }
   }

   // Perform log-rank test
   double hr = calculateHazardRatio(survivalTimes, enrollmentTimes, exitTimes, AtEventCount{ eventsFaAll });
   if(hr > 1.0)
      return {true, 1.0};
   auto lrStatsFa = computeLogRankComponents(survivalTimes, enrollmentTimes, exitTimes, AtEventCount{ eventsFaAll });
   auto chiSquare = logRankStatistic(lrStatsFa);
   if (lrStatsFa.totalEvents != eventsFaAll) {
      std::cout << "Error, events Fa: " << lrStatsFa.totalEvents << std::endl;
      abort();
   }
   double pValue = chiSquaredToPvalue(chiSquare, 1) / 2.0;

   return {true, pValue};
}

int main() {
   constexpr std::array<size_t, arms> numTestsControl = { 25, 25 };
   auto calcMeadianSoBat = [](size_t i) { return i * 1.0 + 5.0; };
   auto calcMeadianSoGps = [](size_t i) { return i * 1.0 + 5.0; };

   // shape Raster: 0.4 .. 1.0 in 0.1 steps
   constexpr auto createShapeData = [](){
      if constexpr (shapeGrid) {
         constexpr size_t numShape = 7;
         return std::array<double, numShape>{0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
      }
      else {
         constexpr size_t numShape = 1;
         return std::array<double, numShape>{shapeDefault};
      }
   };
   auto shapeVals = createShapeData();
   constexpr size_t numShape = shapeVals.size();

   // Calculate mu and sigma for log-normal distribution
   double mu_control, sigma_control, mu_treatment, sigma_treatment;
   lognormal_params(median_control_prior, std_control_prior, &mu_control, &sigma_control);
   lognormal_params(median_treatment_prior, std_treatment_prior, &mu_treatment, &sigma_treatment);

   // Counter for successful studies and weighted sum (4D grids: medianBat x medianGps x shapeBat x shapeGps)
   std::array<std::array<std::array<std::array<double, numShape>, numShape>, numTestsControl[Arm::Bat]>, numTestsControl[Arm::Gps]> successRate = {};
   std::array<std::array<std::array<std::array<double, numShape>, numShape>, numTestsControl[Arm::Bat]>, numTestsControl[Arm::Gps]> validRate = {};

   // Calculate weights for each combination (4D prior)
   std::array<std::array<std::array<std::array<double, numShape>, numShape>, numTestsControl[Arm::Bat]>, numTestsControl[Arm::Gps]> pMosPrior{};
   double total_weight = 0.0;
   double prior_shape = 1.0 / (static_cast<double>(numShape) * static_cast<double>(numShape)); // uniform over shape pairs

   for (size_t c = 0; c < numTestsControl[Arm::Bat]; c++) {
      for (size_t t = 0; t < numTestsControl[Arm::Gps]; t++) {
         double prior_control = lognormal_pdf(calcMeadianSoBat(c), mu_control, sigma_control);
         double prior_treatment = lognormal_pdf(calcMeadianSoGps(t), mu_treatment, sigma_treatment);
         for (size_t sb = 0; sb < numShape; ++sb) {
            for (size_t sg = 0; sg < numShape; ++sg) {
               pMosPrior[c][t][sb][sg] = prior_control * prior_treatment * prior_shape;
               total_weight += pMosPrior[c][t][sb][sg];
            }
         }
      }
   }

   if (total_weight <= 0.0) {
      std::cerr << "Error: Total weight is 0, no valid priors." << std::endl;
      return 1;
   }

   for (size_t c = 0; c < numTestsControl[Arm::Bat]; c++) {
      for (size_t t = 0; t < numTestsControl[Arm::Gps]; t++) {
         for (size_t sb = 0; sb < numShape; ++sb) {
            for (size_t sg = 0; sg < numShape; ++sg) {
               pMosPrior[c][t][sb][sg] /= total_weight;
            }
         }
      }
   }

   const std::array<size_t, arms> numPatients = { maxPatientsPerArm, maxPatientsPerArm };
   // Simulation über 4D grid
   for (size_t c = 0; c < numTestsControl[Arm::Bat]; c++) {
      for (size_t t = 0; t < numTestsControl[Arm::Gps]; t++) {
         for (size_t sb = 0; sb < numShape; ++sb) {
            for (size_t sg = 0; sg < numShape; ++sg) {

               std::array<double, arms> medianSurvival;
               medianSurvival[Arm::Bat] = calcMeadianSoBat(c);
               medianSurvival[Arm::Gps] = calcMeadianSoGps(t);

               std::array<double, arms> shape;
               shape[Arm::Bat] = shapeVals[sb];
               shape[Arm::Gps] = shapeVals[sg];

               size_t successCount = 0;
               size_t simulationsCount = 0;
               for (size_t sim = 0; sim < numSimulations; sim++) {
                  auto [valid, pValue] = simulateInstance(
                     medianSurvival,
                     shape,
                     numPatients
                  );

                  if (valid && pValue < pLimitFa)
                     successCount++;
                  if (valid)
                     simulationsCount++;
               }

               successRate[c][t][sb][sg] = (simulationsCount > 0) ? static_cast<double>(successCount) / static_cast<double>(simulationsCount) : 0.0;
               validRate[c][t][sb][sg] = (numSimulations > 0) ? static_cast<double>(simulationsCount) / static_cast<double>(numSimulations) : 0.0;

               std::cout << "progress: " << std::fixed << std::setprecision(4)
                  << static_cast<double>((sg + sb * numShape) + (t + c * numTestsControl[Arm::Gps]) * (numShape * numShape))
                  / ((numTestsControl[Arm::Bat] * numTestsControl[Arm::Gps]) * (numShape * numShape))
                  << std::endl;
            }
         }
      }
   }

   // pValid (marginal probability that the scenario is "valid")
   double pValid = 0.0;
   for (size_t c = 0; c < numTestsControl[Arm::Bat]; c++) {
      for (size_t t = 0; t < numTestsControl[Arm::Gps]; t++) {
         for (size_t sb = 0; sb < numShape; ++sb) {
            for (size_t sg = 0; sg < numShape; ++sg) {
               pValid += validRate[c][t][sb][sg] * pMosPrior[c][t][sb][sg];
            }
         }
      }
   }

   // Posterior on 4D grid
   std::array<std::array<std::array<std::array<double, numShape>, numShape>, numTestsControl[Arm::Bat]>, numTestsControl[Arm::Gps]> pMosPost = {};
   for (size_t c = 0; c < numTestsControl[Arm::Bat]; c++) {
      for (size_t t = 0; t < numTestsControl[Arm::Gps]; t++) {
         for (size_t sb = 0; sb < numShape; ++sb) {
            for (size_t sg = 0; sg < numShape; ++sg) {
               if (pValid > 0.0)
                  pMosPost[c][t][sb][sg] = validRate[c][t][sb][sg] * pMosPrior[c][t][sb][sg] / pValid;
               else
                  pMosPost[c][t][sb][sg] = 0.0;
            }
         }
      }
   }

   std::cout << "Parameters of the log-normal priorities:" << std::endl;
   std::cout << "Control: E[X] = " << std::fixed << std::setprecision(1) << median_control_prior
      << ", Std[X] = " << std_control_prior
      << " -> mu = " << std::fixed << std::setprecision(4) << mu_control
      << ", sigma = " << sigma_control << std::endl;
   std::cout << "Treatment: E[X] = " << std::fixed << std::setprecision(1) << median_treatment_prior
      << ", Std[X] = " << std_treatment_prior
      << " -> mu = " << std::fixed << std::setprecision(4) << mu_treatment
      << ", sigma = " << sigma_treatment << std::endl;

   std::ofstream csvFile("simulation_results.csv");
   if (!csvFile.is_open()) {
      std::cerr << "error while opening result file" << std::endl;
      return 1;
   }

   csvFile << "MedianBat,MedianGps,ShapeBat,ShapeGps,SuccessRate,ValidRate,PosteriorWeight\n";
   for (size_t sb = 0; sb < numShape; ++sb) {
      for (size_t sg = 0; sg < numShape; ++sg) {
         for (size_t c = 0; c < numTestsControl[Arm::Bat]; ++c) {
            for (size_t t = 0; t < numTestsControl[Arm::Gps]; ++t) {
               double medianBat = calcMeadianSoBat(c);
               double medianGps = calcMeadianSoGps(t);

               csvFile
                  << medianBat << ","
                  << medianGps << ","
                  << shapeVals[sb] << ","
                  << shapeVals[sg] << ","
                  << successRate[c][t][sb][sg] << ","
                  << validRate[c][t][sb][sg] << ","
                  << pMosPost[c][t][sb][sg]
                  << "\n";
            }
         }
      }
   }

   // Total weighted propability of success
   double pSuccess = 0.0;
   double pMosPostSum = 0.0;
   for (size_t c = 0; c < numTestsControl[Arm::Bat]; c++) {
      for (size_t t = 0; t < numTestsControl[Arm::Gps]; t++) {
         for (size_t sb = 0; sb < numShape; ++sb) {
            for (size_t sg = 0; sg < numShape; ++sg) {
               pSuccess += successRate[c][t][sb][sg] * pMosPost[c][t][sb][sg];
               pMosPostSum += pMosPost[c][t][sb][sg];
            }
         }
      }
   }
   pSuccess = (pMosPostSum > 0) ? pSuccess / pMosPostSum : 0.0;
   printf("\nTotal weighted propability of success: %.6f\n", pSuccess);

   return 0;
}
