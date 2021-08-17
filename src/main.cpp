#include "getFoodData.h"
#include "genetic_algo.h"
#include <memory>
#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>


void interpretRule(const std::vector<std::string> &rule,const FoodDatabase &foodDatabase, std::vector<int> &InterpretedRule);
void debugMaps(const std::map<std::string, int> &foodMap);
std::vector<float> getDecisionVariableConstrain();
int main()
{
    std::string FileLocation = "/home/abhi/ga-cpp/files/food_data_2.csv";
    getFoodData food_data(FileLocation);
    FoodDatabase foodDB = food_data.getFoodDataList(); 
    std::vector<std::string> rule = {"Main","Veggies","Veggies","Side"};  //.getRule()
    std::vector<int> InterpretedRule;
    InterpretedRule.clear();
    InterpretedRule.reserve(rule.size());
    interpretRule(rule,foodDB,InterpretedRule);
    std::vector<float> decisionVariableConstrain = getDecisionVariableConstrain();

    int populationSize = 1000; 
    float idealCalorie = 700; //a return value of a customer class (Struct customer-cal, prot, fat, carbs, liking, dislike )
    uint numberOfIteration = 200;
    uint8_t numberOfRecommendedMeal = 2;
    float servingSizeMultipler = 0.1f;

    std::unique_ptr<geneticAlgo> genAlgo (new geneticAlgo());
    genAlgo->setPopulationSize(populationSize);
    genAlgo->setIdealCalorie(idealCalorie);
    genAlgo->setRuleComposition(InterpretedRule);
    genAlgo->setFoodDatabase(foodDB);
    genAlgo->setNSGAbasedOptimization(true);
    genAlgo->setNumGeneration(numberOfIteration);
    genAlgo->setServingSizeParts(servingSizeMultipler);
    std::vector<geneticAlgo::Meal> recommenedMeal = genAlgo->getRecommendedMeal(numberOfRecommendedMeal);
    std::cout << recommenedMeal.size() << std::endl;

    return 0;
}

std::vector<float> getDecisionVariableConstrain()
{
    std::vector<float> decisionVar;
    decisionVar.reserve(4);
    decisionVar[inverseCalDiff] = 1.0F; //
    decisionVar[carbCal] = 1.0F; // 
    decisionVar[proteinCal] = 1.0F; //
    decisionVar[FatCal] = 1.0F; //
    return decisionVar;
}

void interpretRule(const std::vector<std::string> &rule,const FoodDatabase &foodDatabase,std::vector<int> &InterpretedRule)
{
    for (const std::string &ruleComp : rule)
    {
        InterpretedRule.push_back(foodDatabase.foodCategories.at(ruleComp));
    }
}

void debugMaps(const std::map<std::string, int> &foodMap)
{
    for (auto const&x : foodMap)
    {
        std::cout << x.first << ":" << x.second <<std::endl;
    }
}