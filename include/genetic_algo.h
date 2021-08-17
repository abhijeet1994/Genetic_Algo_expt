#ifndef GENETIC_ALGO_H
#define GENETIC_ALGO_H
#include "getFoodData.h"
#include <array>
#include <limits>

namespace geneticAlgoConsts
{
    constexpr float MultiConstProtein = 4.0f; 
    constexpr float MinCalProtein = 0.10f;
    constexpr float MaxCalProtein = 0.35f;
    constexpr float MultiConstCarb = 4.0f;
    constexpr float MinCalCarb = 0.45f;
    constexpr float MaxCalCarb = 0.65f;
    constexpr float MultiConstFat = 9.0F;
    constexpr float MinCalFat = 0.20f;
    constexpr float MaxCalFat = 0.35f;
    constexpr uint  DefaultPopulationSize = 1000;
    constexpr uint DefaultNumIteration = 75;
    constexpr float ServingSizeParts = 0.1f;
}

enum decisionVariables {inverseCalDiff=0, proteinCal, carbCal, FatCal};

class geneticAlgo
{
    public:
    // geneticAlgo(const int &population, const FoodDatabase &foodDatabase,const std::vector<int> &interpretedRule,const float &idealCalorie,const std::vector<float> &decisionVariableConstrain);
    geneticAlgo();
    void recommendMeal();
    void recommendMealUsingNsga();
    struct Meal
    {
        std::vector<int> mealComponent;
        std::vector<float> servingSize;
        float calorieInMeal;
        float protein_cont;
        float carbs_cont;
        float fats_cont;
        float fitness_value;
        std::array<float,2> fitness_cpd; //cumulative probability distribution
        std::vector<float> fitness_array;
    };

    struct MealFitness 
    {
        std::vector<float> fitnessArray;
        int MealIndex;
    };

    struct MealDistances
    {
        Meal meal;
        int ind_index;
        float distance;
    };

    //setter function
    void setPopulationSize(const int &populationSize);
    void setIdealCalorie(const int &idealCalorie);
    void setRuleComposition(const std::vector<int> &interpretedRule);
    void setFoodDatabase(const FoodDatabase &foodDatabase);
    void setNSGAbasedOptimization(const bool &isNSGA);
    void setCarbCalorieProportion(const float &minCarbCal, const float &maxCarbCal);
    void setProteinCalorieProportion(const float &minProteinCal, const float &maxProteinCal);
    void setFatCalorieProportion(const float &minFatCal, const float &maxFatCal);
    void setNumGeneration(const uint &numIteration);
    void setServingSizeParts(const float &servingSizeMultiper);
    std::vector<Meal> getRecommendedMeal(const uint &numberOfMeals);
    //
    private:
    void createMealPopulationStruct();
    void fitnessBasedMealStructSort(std::vector<Meal> &mealVector);
    void selectFirstNMeal(std::vector<Meal> &mealVector, const int &NumberToSelect);
    std::vector<int> generateRandomIntVector(const int &lastIndex, const int &sizeOfVector );
    float generateRandomServingSize();
    int binaryProbSearch(const std::vector<Meal> &mealPopulation,const float &probability);
    int generateRandomInt(const int &lastIndex);
    float generate5thDigitProb();
    float createCompleteMealCrossoverPopulation(std::vector<Meal> &mealVector);
    void createMealCrossoverPopulation(std::vector<Meal> &mealVector);
    
    void calFitnessProbForEachMeal(std::vector<Meal> &mealVector,const float &fitness_sum);
    void mutation(Meal &child);
    void mutateMealComponent(std::vector<int> &mealComp);
    void mutateServingSize(std::vector<float> &servingSize);
    void calculateMealMacro(Meal &meal);
    void calculateSelectionProbForMealVec(std::vector<Meal> &mealVector,const float &fitness_sum);
    void combineMealVectors(std::vector<Meal> &Meal_1, std::vector<Meal> &Meal_2, std::vector<Meal> &outputMeal);
    void copyMealVector(const std::vector<Meal> &originalVector, std::vector<Meal> &copyVector);
    float calculateFitnessSum(std::vector<Meal> &mealVector);

    float getMacroFitness (const float &givenMacro, const float &minMacroReq, const float &maxMacroAllowed);
    void updateMealFitnessArray( Meal &meal);
    void sortNonDominated();
    int nsgaCompareFast(const std::vector<float> &a, const std::vector<float> &b, const std::vector<float> &sign);
    void calculateMealDistance();
    void sortNsgaPopulation();
    void initializingFitnessArray();
    void updateFitnessArray(const uint &sizeFitnessArray, const uint &sizeFitnessSign);
    void classifyMealDomination();
    void createParetoFronts();
    void initializeMealFitnessPopulation();
    void initialzeMealDistancePopulation(const std::vector<int> &front);
    
    uint _mPopulationSize;
    float _mFitnessSum;
    int _mcounter;
    bool _mIsNSGA;
    uint _mCombinedPopulationSize;
    std::vector<float> _mFitness_sign;
    float _mMinCalProtein;
    float _mMaxCalProtein;
    float _mMinCalCarb;
    float _mMaxCalCarb;
    float _mMinCalFat;
    float _mMaxCalFat;
    uint _mNumCompoInRule;
    uint _mNumberIteration;
    float _mServingSizeParts;
    std::vector<int> _mCompPositionInRule;
    std::vector<std::vector<int>> _mMealItems;
    std::vector<std::vector<float>> _mServingSize;
    std::vector<std::vector<foodData>> _mfoodList;
    std::vector<float> _mCalorieInMeal;
    std::vector<std::vector<float>> _mMacrosInMeal;

    std::vector<Meal> _mMealPopulation;
    std::vector<Meal> _mMealCrossoverPopulation;
    std::vector<Meal> _mCombinedMealPopulation;
    std::vector<Meal> _mTempMealPopulation;
    std::vector<Meal> _mCopyTempMeal;
    std::map<std::string, int> _mFoodCategories;
    std::vector<std::vector<float>> _mfitnessArray;
    std::vector<std::vector<int>> _mfronts;
    std::vector<float> _mDistances;

    std::vector<int> _mCurrentFront;
    std::vector<int> _mNextFront;
    std::vector<int> _mDominatingMeals; //Stores number of meals dominating a given ith meal
    std::vector<std::vector<int>> _mDominatedMeals; //Meals which ith meal dominate
    std::vector<MealFitness> _mMealFitnessPopulation;
    std::vector<MealDistances> _mNextGroupPopulation;
    float _idealCalorie;
};

#endif 