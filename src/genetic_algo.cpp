#include "genetic_algo.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <random>

geneticAlgo::geneticAlgo():
_mPopulationSize(geneticAlgoConsts::DefaultPopulationSize),
_mFitnessSum(0.0),
_mcounter(0),
_mIsNSGA(false),
_mCombinedPopulationSize(0),
//_mFitness_sign({1.0F,1.0F,1.0F,1.0F}),
_mFitness_sign({1.0F}),
_mMinCalProtein(geneticAlgoConsts::MinCalProtein),
_mMaxCalProtein(geneticAlgoConsts::MaxCalProtein),
_mMinCalCarb(geneticAlgoConsts::MinCalCarb),
_mMaxCalCarb(geneticAlgoConsts::MaxCalCarb),
_mMinCalFat(geneticAlgoConsts::MinCalFat),
_mMaxCalFat(geneticAlgoConsts::MaxCalFat),
_mNumberIteration(geneticAlgoConsts::DefaultNumIteration),
_mServingSizeParts(geneticAlgoConsts::ServingSizeParts)
{
    _mMealPopulation.reserve(_mPopulationSize);
    _mTempMealPopulation.reserve(_mPopulationSize+_mPopulationSize);
    _mCopyTempMeal.reserve(_mPopulationSize+_mPopulationSize);
    _mMealCrossoverPopulation.reserve(_mPopulationSize);
    _mfitnessArray.reserve(_mPopulationSize+_mPopulationSize);
    _mCurrentFront.reserve(_mPopulationSize);
    _mNextFront.reserve(_mPopulationSize);
    _mDominatingMeals.reserve(_mPopulationSize+_mPopulationSize);
    _mDominatedMeals.reserve(_mPopulationSize+_mPopulationSize);
    _mMealFitnessPopulation.reserve(_mPopulationSize+_mPopulationSize);
    _mDistances.reserve(_mPopulationSize+_mPopulationSize);
    _mNextGroupPopulation.reserve(_mPopulationSize+_mPopulationSize);
}
void geneticAlgo::setPopulationSize(const int &populationSize)
{
    _mPopulationSize = populationSize;
}

void geneticAlgo::setIdealCalorie(const int &idealCalorie)
{
    _idealCalorie = idealCalorie;
}
void geneticAlgo::setRuleComposition(const std::vector<int> &interpretedRule)
{
    _mCompPositionInRule = interpretedRule;
    _mNumCompoInRule = interpretedRule.size();
}
void geneticAlgo::setFoodDatabase(const FoodDatabase &foodDatabase)
{
    _mfoodList = foodDatabase.food_database;
    _mFoodCategories = foodDatabase.foodCategories;
}
void geneticAlgo::setNSGAbasedOptimization(const bool &isNSGA)
{
    _mIsNSGA = isNSGA;
}
void geneticAlgo::setCarbCalorieProportion(const float &minCarbCal, const float &maxCarbCal)
{
    _mMinCalCarb=minCarbCal;
    _mMaxCalCarb=maxCarbCal;
}
void geneticAlgo::setProteinCalorieProportion(const float &minProteinCal, const float &maxProteinCal)
{
    _mMinCalProtein=minProteinCal;
    _mMaxCalProtein=maxProteinCal;
}
void geneticAlgo::setFatCalorieProportion(const float &minFatCal, const float &maxFatCal)
{
    _mMinCalFat=minFatCal;
    _mMaxCalFat=maxFatCal;
}
void geneticAlgo::setNumGeneration(const uint &numIteration)
{
    _mNumberIteration = numIteration;
}
void geneticAlgo::setServingSizeParts(const float &servingSizeMultiper)
{
    _mServingSizeParts = servingSizeMultiper;
}

std::vector<geneticAlgo::Meal> geneticAlgo::getRecommendedMeal(const uint &numberOfMeals)
{
    if (_mIsNSGA == true)
    {
        recommendMealUsingNsga();
    }
    else
    {
        recommendMeal();
    }
    _mMealPopulation.resize(numberOfMeals);
    return _mMealPopulation;

}
void geneticAlgo::recommendMealUsingNsga()
{
    _mIsNSGA = true;
    float fitness_sum;
    createMealPopulationStruct();
    // uint numberIteration=50;
    for (size_t iterator=0; iterator < _mNumberIteration; ++iterator)
    {
        fitness_sum =   createCompleteMealCrossoverPopulation(_mMealPopulation);   
        _mTempMealPopulation.clear();
        combineMealVectors(_mMealPopulation, _mMealCrossoverPopulation,_mTempMealPopulation);
        sortNonDominated();
        calculateMealDistance();
        sortNsgaPopulation();
        selectFirstNMeal(_mTempMealPopulation, _mPopulationSize);
        _mMealPopulation.clear();
        copyMealVector(_mTempMealPopulation, _mMealPopulation);
        

        if (iterator % 20 == 0)
        {
            std::cout << _mMealPopulation[1].calorieInMeal <<", FV ," <<_mMealPopulation[1].fitness_value<<" , iteration : "<<iterator <<std::endl;
            for (int i = 0; i < 10; i++)
            {
                std::cout<< "individual rank : " << i << " : ";
                std::cout << _mMealPopulation[i].calorieInMeal << " , ";
                for (int j = 0 ; j < _mMealPopulation[0].fitness_array.size(); j ++)
                {
                    std::cout << _mMealPopulation[i].fitness_array[j] << " , ";
                }
                std::cout << "\n" ;
            }
        }
    }

    for (size_t i=0; i < 2; i++)
    {
    std::vector<int> bestMealComponent =  _mMealPopulation[i*5].mealComponent;
    std::vector<float> bestServingSize =  _mMealPopulation[i*5].servingSize;  
    std::cout << "Food Rank: " << i*5 << std::endl;
    std::cout << "caloire:" << _mMealPopulation[i*5].calorieInMeal << " protein ratio:" << (_mMealPopulation[i*5].protein_cont * 4.0) / _mMealPopulation[i*5].calorieInMeal<< " carb ratio:" << (_mMealPopulation[i*5].carbs_cont * 4.0) / _mMealPopulation[i*5].calorieInMeal << std::endl;
    for (uint j=0; j<_mNumCompoInRule; ++j)
    {
        std::cout <<"FOOD NAME: " << _mfoodList.at(_mCompPositionInRule[j])[bestMealComponent[j]].ingredientName << ", SS: " << bestServingSize[j] <<std::endl;
    }
        std::cout << "----------------------------" << std::endl;
    }
}

void geneticAlgo::sortNsgaPopulation()
{
    _mCopyTempMeal.clear();
    _mCopyTempMeal = _mTempMealPopulation;
    _mTempMealPopulation.clear();
    _mNextGroupPopulation.clear();

    // std::cout << "Befor sorted fronts : ";
    // for (std::vector<int> &front : _mfronts)
    // {
    //     for (int &ind : front)
    //     {
    //         std::cout << _mCopyTempMeal[ind].calorieInMeal<< " , ";
    //     }         
    // }
    // std::cout<< "\n";
    
    for (std::vector<int> &front : _mfronts)
    {
        initialzeMealDistancePopulation(front);
        std::sort(_mNextGroupPopulation.begin(), _mNextGroupPopulation.end(),[&](const MealDistances & a, const MealDistances & b)
        {
            return a.distance > b.distance;
        });
        for (MealDistances &nextMeal : _mNextGroupPopulation)
        { 
            _mTempMealPopulation.push_back(nextMeal.meal);
        }
    }
    _mCopyTempMeal.clear();
}

void geneticAlgo::initialzeMealDistancePopulation(const std::vector<int> &front)
{
    for (const int individual : front)
        {
            MealDistances tmp;
            tmp.meal = _mCopyTempMeal[individual];
            tmp.ind_index = individual;
            tmp.distance = _mDistances[individual];
            _mNextGroupPopulation.push_back(tmp);
        }
}

void geneticAlgo::sortNonDominated()
{
    bool first_front_only = false ; // Only sorts the first front. We want to sort all individuals. Thus we set it false
    // std::vector<float> _mFitness_sign{1.0F, 1.0F, 1.0F, 1.0F}; // SUPER IMPORTANT -> INPUT FROM OUTSIDE

    _mCombinedPopulationSize = _mTempMealPopulation.size();
    _mfitnessArray.clear();
    _mCurrentFront.clear(); //Population indices which are in current front
    _mNextFront.clear();
    _mDominatingMeals.clear(); //No of individuals dominating you//
    _mDominatedMeals.clear(); ////Index of people you dominate////
    //Population[x] dominates Population[y],Population[a],Population[b],Population[c]
    //dominated[x].get[] -->Ans [a,b,c] 

    initializingFitnessArray();
    updateFitnessArray(_mfitnessArray.size(),_mFitness_sign.size());
    classifyMealDomination(); // Creating a list of no of individuals dominating and the individuals you dominate 
    
    if (first_front_only == false)
    {
        createParetoFronts();
    }
    /* Debug
    int chk =0;
    for (std::vector<int> &front : _mfronts){ 
        for (int &indi : front){
            chk += 1;}}
    if (chk < _mCombinedPopulationSize){
        std::cout << " ******** Error Detected At Check ********* " << std::endl;
        std::cout << "Current_front : " << _mCurrentFront.size() <<std::endl;
        std::cout << chk << std::endl;
    }   
    */
}

void geneticAlgo::initializingFitnessArray()
{
    std::vector<int> tmp;
    tmp.clear();
    for (size_t iterator = 0; iterator < _mCombinedPopulationSize; ++iterator )
    {
        _mfitnessArray.push_back(_mTempMealPopulation[iterator].fitness_array);
        _mDominatingMeals.push_back(0);
        _mDominatedMeals.push_back(tmp);
    }
}

void geneticAlgo::updateFitnessArray(const uint &sizeFitnessArray, const uint &sizeFitnessSign)
{
    for (size_t j=0; j < sizeFitnessArray; ++j)
    {
        for (size_t i=0; i <sizeFitnessSign; ++i)
        {
            _mfitnessArray[j][i] = _mfitnessArray[j][i]*_mFitness_sign[i]; 
        }
    }
}

void geneticAlgo::classifyMealDomination()
{
    for (size_t i = 0 ; i < _mCombinedPopulationSize ; ++i)
    {
        for (size_t j = i + 1 ; j < _mCombinedPopulationSize; ++j)
        {    
            int tmp = nsgaCompareFast(_mfitnessArray[i], _mfitnessArray[j], _mFitness_sign);
            if (tmp == 1)
            {
                _mDominatingMeals[j] += 1;
                _mDominatedMeals[i].push_back(j);
            }
            else if (tmp == -1)
            {
                _mDominatingMeals[i] += 1;
                _mDominatedMeals[j].push_back(i);
            }
        }
        if ( _mDominatingMeals[i] == 0)
        {
            _mCurrentFront.push_back(i);
            //Meal which are not dominated by any meals as of this instant gets added to current front
        }
    } //List of dominating fits, dominated fits, and current_front created.
}

void geneticAlgo::createParetoFronts()
{
    _mfronts.clear();
    uint pareto_sorted = _mCurrentFront.size(); //Individuals sorted
    while (pareto_sorted < _mCombinedPopulationSize || _mCurrentFront.size() != 0 )
    {
        std::vector<int> tmp; 
        _mfronts.push_back(tmp);
        //_mfront[0] --> holds indices of all 1st best meal
        for (int &individual : _mCurrentFront)
        {
            for (int &dom_ind : _mDominatedMeals[individual])
            {
                    _mDominatingMeals[dom_ind] -= 1; 
                if ( _mDominatingMeals[dom_ind] ==0 ) //Meals which are not being dominated than added to next front (2nd best meals)
                {
                    _mNextFront.push_back(dom_ind); //Individual gettings sorted (next fornt)
                    pareto_sorted += 1; //individual sorted gets incremented
                }  
            }
        }
        _mfronts[_mfronts.size() - 1] = _mCurrentFront;
        _mCurrentFront = _mNextFront;
        /* Debug
        if (_mCurrentFront.size() == 0 && pareto_sorted < _mCombinedPopulationSize)
        {
            std::cout << " ******** Error Detected ********* " << std::endl;
        }
        */
        _mNextFront.clear();
    }
}

void geneticAlgo::calculateMealDistance()
{
    _mDistances.clear();
    _mMealFitnessPopulation.clear();
    initializeMealFitnessPopulation();
    uint fitness_size = _mfitnessArray[0].size();

    for (size_t iterator =0; iterator < fitness_size; ++iterator )
    {
        std::sort(_mMealFitnessPopulation.begin(), _mMealFitnessPopulation.end(), [&](const MealFitness &a, const MealFitness &b)
        {
            return a.fitnessArray[iterator] < b.fitnessArray[iterator];
        });
    
    
    _mDistances[_mMealFitnessPopulation[0].MealIndex] = std::numeric_limits<float>::infinity();
    _mDistances[_mMealFitnessPopulation[_mCombinedPopulationSize-1].MealIndex] = std::numeric_limits<float>::infinity();

    if (_mMealFitnessPopulation[0].fitnessArray[iterator] != _mMealFitnessPopulation[_mCombinedPopulationSize-1].fitnessArray[iterator])
        {
            float norm = _mMealFitnessPopulation[_mCombinedPopulationSize - 1].fitnessArray[iterator] - _mMealFitnessPopulation[_mCombinedPopulationSize - 1].fitnessArray[iterator];
            for (size_t j=1; j < _mCombinedPopulationSize-1; ++j)
            {
                int curr_ind = _mMealFitnessPopulation[j].MealIndex;
                float dist_smaller = _mMealFitnessPopulation[j - 1].fitnessArray[iterator];
                float dist_greater = _mMealFitnessPopulation[j + 1].fitnessArray[iterator];
                _mDistances[curr_ind] += (dist_greater - dist_smaller)/norm;
            }
        }
    }
}

void geneticAlgo::initializeMealFitnessPopulation()
{

    for (size_t iterator =0; iterator < _mCombinedPopulationSize; ++iterator)
    {
        MealFitness ind;
        ind.fitnessArray = _mfitnessArray[iterator];
        ind.MealIndex = iterator;
        _mMealFitnessPopulation.push_back(ind);
    }

}

int geneticAlgo::nsgaCompareFast(const std::vector<float> &a, const std::vector<float> &b, const std::vector<float> &sign)
{
    //if b is better than a in all 4 characteristic than return -1
    //If a is better than b in all 4 characteristic than return +1
    //If none of the above holds true pass 0
    bool flag = false;
    // int dominate = 0;
    uint32_t sizeSign = sign.size();
    for (size_t iterator =0; iterator < sizeSign; ++iterator)
    {
        float tmp = a[iterator] - b[iterator];
        if (tmp ==0) return 0;
        if (iterator == 0) 
        {
            if (tmp > 0 ) flag = true;
            else flag =false;
        }
        if ((tmp >0) !=flag ) return 0;
    }
    if (flag == true) return 1;
    return -1;
}

/**
 * @brief This function gets called for calculating recommedned meal based on closer to ideal calorie
 * 
 */
void geneticAlgo::recommendMeal()
{
    float fitness_sum =0.0;
    createMealPopulationStruct();
    // uint numberIteration = 50;
    for (size_t iterator =0; iterator < _mNumberIteration; ++iterator )
    {
    calculateSelectionProbForMealVec(_mMealPopulation,_mFitnessSum);
    fitness_sum =   createCompleteMealCrossoverPopulation(_mMealPopulation);
    calculateSelectionProbForMealVec(_mMealCrossoverPopulation,fitness_sum);
    _mTempMealPopulation.clear();
    combineMealVectors(_mMealPopulation, _mMealCrossoverPopulation,_mTempMealPopulation);
    fitnessBasedMealStructSort(_mTempMealPopulation);
    selectFirstNMeal(_mTempMealPopulation,_mPopulationSize);
    _mMealPopulation.clear();
    copyMealVector(_mTempMealPopulation,_mMealPopulation);
    fitness_sum = calculateFitnessSum(_mMealPopulation);
    }

    for (size_t i=0; i < 2; i++)
    {
    std::vector<int> bestMealComponent =  _mMealPopulation[i*5].mealComponent;
    std::vector<float> bestServingSize =  _mMealPopulation[i*5].servingSize;  
    std::cout << "Food Rank: " << i*5 << std::endl;
    std::cout << "caloire:" << _mMealPopulation[i*5].calorieInMeal << " protein ratio:" << (_mMealPopulation[i*5].protein_cont * 4.0) / _mMealPopulation[i*5].calorieInMeal<< " carb ratio:" << (_mMealPopulation[i*5].carbs_cont * 4.0) / _mMealPopulation[i*5].calorieInMeal << std::endl;
    for (size_t j=0; j<_mNumCompoInRule; ++j)
    {
        std::cout <<"FOOD NAME: " << _mfoodList.at(_mCompPositionInRule[j])[bestMealComponent[j]].ingredientName << ", SS: " << bestServingSize[j] <<std::endl;
    }
        std::cout << "----------------------------" << std::endl;
    }

}
/**
 * @brief Create initial meal population based on rule from customer. 
 * 
 * @param idealCalorie Ideal required calorie of customer
 */

void geneticAlgo::createMealPopulationStruct()
{
    _mMealPopulation.clear();
    for (size_t i=0; i < _mPopulationSize; ++i)
    {
        Meal meal;
        for (size_t j=0; j < _mNumCompoInRule; ++j)
        {
            meal.mealComponent.push_back(generateRandomInt(_mfoodList[_mCompPositionInRule[j]].size()-1));
            meal.servingSize.push_back(generateRandomServingSize());
        }
        calculateMealMacro(meal);
        _mMealPopulation.push_back(meal);
    }
}
/**
 * @brief Creates crossover population and return fitness sum parameter calculated from crossover population
 * 
 * @param mealVector Parent population from which crossover is generated
 * @return float fitness_sum parameter
 */
float geneticAlgo::createCompleteMealCrossoverPopulation(std::vector<Meal> &mealVector)
{
    _mMealCrossoverPopulation.clear();
    
    while (_mMealCrossoverPopulation.size() < _mPopulationSize)
    {   
        createMealCrossoverPopulation(mealVector);
    }
    //std::cout<<"\n";
    return _mFitnessSum;
}

/**
 * @brief Thsi function generates two crossover child
 * 
 * @param mealVector Parent population 
 */
void geneticAlgo::createMealCrossoverPopulation(std::vector<Meal> &mealVector)
{  
    //std::cout<<"Crashes there "; 
    int parent_1, parent_2;
     
    if (_mIsNSGA == false)
    {
        float parent1 = generate5thDigitProb();
        float parent2 = generate5thDigitProb();
        parent_1 = binaryProbSearch(mealVector,parent1);
        parent_2 = binaryProbSearch(mealVector,parent2);
    }
    else
    {
        parent_1 = std::min(generateRandomInt(_mPopulationSize -1), generateRandomInt(_mPopulationSize -1));
        parent_2 = std::min(generateRandomInt(_mPopulationSize -1), generateRandomInt(_mPopulationSize -1));
    }
      
    //Randomly select how many meal + its corresponding service size to interchange
    Meal child_1 = mealVector[parent_1];
    Meal child_2 = mealVector[parent_2];
    
    //Select number of genes to switch between 2 child [from 0 to size of mealComponent]
    uint geneToswap = generateRandomInt((child_1.mealComponent.size() -1));
    //Select which gene to swap based on number of geneToswap
    std::vector<int> indexToSwap = generateRandomIntVector((child_1.mealComponent.size() -2),geneToswap);
    for (size_t iterator=0; iterator < geneToswap; ++iterator)
    {
        int tempMealComp;
        float tempServingSize;
        tempMealComp = child_1.mealComponent.at(indexToSwap.at(iterator));
        tempServingSize = child_1.servingSize.at(indexToSwap.at(iterator));
        child_1.mealComponent.at(indexToSwap.at(iterator)) = child_2.mealComponent.at(indexToSwap.at(iterator));
        child_1.servingSize.at(indexToSwap.at(iterator)) = child_2.servingSize.at(indexToSwap.at(iterator));
        child_2.mealComponent.at(indexToSwap.at(iterator)) = tempMealComp;
        child_2.servingSize.at(indexToSwap.at(iterator)) = tempServingSize;
    }   

    mutation(child_1);
    mutation(child_2);
    
    calculateMealMacro(child_1);
    calculateMealMacro(child_2);
    
    _mMealCrossoverPopulation.push_back(child_1); 
    _mMealCrossoverPopulation.push_back(child_2);
    //_mcounter++;

}
/**
 * @brief Combines two meal population 
 * 
 * @param Meal_1 Meal population 1
 * @param Meal_2 Meal population 2
 * @param outputMeal Combined meal population
 */
void geneticAlgo::combineMealVectors(std::vector<Meal> &Meal_1, std::vector<Meal> &Meal_2, std::vector<Meal> &outputMeal)
{
    outputMeal.clear();
    outputMeal.reserve(Meal_1.size() + Meal_2.size());
    outputMeal.insert(outputMeal.end(),Meal_1.begin(), Meal_1.end());
    outputMeal.insert(outputMeal.end(),Meal_2.begin(), Meal_2.end());
}
/**
 * @brief 
 * 
 * @param mealVector 
 * @param fitness_sum 
 */
void geneticAlgo::calculateSelectionProbForMealVec(std::vector<Meal> &mealVector,const float &fitness_sum)
{
    fitnessBasedMealStructSort(mealVector);
    calFitnessProbForEachMeal(mealVector,fitness_sum);
}
/**
 * @brief Creates copy of vector
 * 
 * @param originalVector (>_<)
 * @param copyVector (> = <)
 */
void geneticAlgo::copyMealVector(const std::vector<Meal> &originalVector, std::vector<Meal> &copyVector)
{
    std::copy(originalVector.begin(),originalVector.end(), std::back_inserter(copyVector));
}
/**
 * @brief Performs mutation on child meal
 * 
 * @param child Meal who's going to be mutate
 */
//int 0->No mutation, 1-> Mutate Meal Component, 2-> mutate serving size
void geneticAlgo::mutation(Meal &child)
{
    int whatTomutate = generateRandomInt(2);
    if (whatTomutate == 0)
    {
        //Do nothing
    }
    else if (whatTomutate ==1)
    {
        mutateMealComponent(child.mealComponent);
    }
    else if (whatTomutate ==2)
    {
        mutateServingSize(child.servingSize);
    }
}
/**
 * @brief Mutate meal component vector
 * 
 * @param mealComp 
 */
void geneticAlgo::mutateMealComponent(std::vector<int> &mealComp)
{
    int whereToMutate = generateRandomInt(mealComp.size() -1);
    int newFoodComp = generateRandomInt(_mfoodList[_mCompPositionInRule[whereToMutate]].size()-1);
    mealComp.at(whereToMutate) = newFoodComp;
}
/**
 * @brief Perform mutation on serving size
 * 
 * @param servingSize Given serving size vector
 */
void geneticAlgo::mutateServingSize(std::vector<float> &servingSize)
{
    int whereToMutate = generateRandomInt(servingSize.size() -1);
    float newServingSize = generateRandomServingSize();
    servingSize.at(whereToMutate) = newServingSize;
}
/**
 * @brief Calculate macros for an given input meal structure 
 * 
 * @param meal Meal structure who macros are to be calculated
 */
void geneticAlgo::calculateMealMacro(Meal &meal)
{
    float fitness_stage1=0.0;
    float calorie = 0.0;
    float protein =0.0;
    float carbs = 0.0;
    float fats = 0.0;
    for (size_t j=0;j<_mNumCompoInRule; ++j)
        {
            calorie = calorie + ((_mfoodList[_mCompPositionInRule[j]][meal.mealComponent[j]].kiloCalories)*meal.servingSize[j]);
            protein = protein + ((_mfoodList[_mCompPositionInRule[j]][meal.mealComponent[j]].proteinContent)*meal.servingSize[j]);
            carbs = carbs + ((_mfoodList[_mCompPositionInRule[j]][meal.mealComponent[j]].carbsContent)*meal.servingSize[j]);
            fats = fats + ((_mfoodList[_mCompPositionInRule[j]][meal.mealComponent[j]].FatContent)*meal.servingSize[j]);
        }
    meal.calorieInMeal = calorie;
    meal.protein_cont = protein;
    meal.carbs_cont = carbs;
    meal.fats_cont = fats;
    fitness_stage1 = (meal.calorieInMeal - _idealCalorie)*(meal.calorieInMeal - _idealCalorie);
    if (fitness_stage1 < 1)
    {
        fitness_stage1 = 1;
    }
    meal.fitness_value = 1000000 * 1/fitness_stage1;
    if (_mIsNSGA == true)
    {
        updateMealFitnessArray(meal);
    }
    _mFitnessSum = _mFitnessSum + meal.fitness_value;
}

void geneticAlgo::updateMealFitnessArray( Meal &meal)
{
    meal.fitness_array.clear();
    meal.fitness_array.push_back(meal.fitness_value*meal.fitness_value);
    
    float proteinCal = (meal.protein_cont * geneticAlgoConsts::MultiConstProtein)/meal.calorieInMeal;
    float carbCal = (meal.carbs_cont * geneticAlgoConsts::MultiConstCarb)/meal.calorieInMeal;
    //float fatCal = (meal.fats_cont * geneticAlgoConsts::MultiConstFat)/meal.calorieInMeal;
    
    meal.fitness_array.push_back(getMacroFitness(proteinCal,geneticAlgoConsts::MinCalProtein,geneticAlgoConsts::MaxCalProtein));
    meal.fitness_array.push_back(getMacroFitness(carbCal,geneticAlgoConsts::MinCalCarb,geneticAlgoConsts::MaxCalCarb));
   // meal.fitness_array.push_back(getMacroFitness(fatCal,geneticAlgoConsts::MinCalFat,geneticAlgoConsts::MaxCalFat));
    
}

float geneticAlgo::getMacroFitness (const float &givenMacro, const float &minMacroReq, const float &maxMacroAllowed)
{    
    //if ((givenMacro >= minMacroReq) && (givenMacro <= maxMacroAllowed)) { return 1.0F;}
    float mean = (maxMacroAllowed + minMacroReq)/2.0 ;
    float val = 1 - 10*(givenMacro - mean)*(givenMacro - mean);
    return val;
}

/**
 * @brief Calculate fitness_sum of a population
 * 
 * @param mealVector Population whose fitness sum is to be calculated
 * @return float Fitness_sum of a population. Consider it as total sum of good qualities of population. The indivudal meal which has 
 * good prperties will have high value of fitness and hence high value of fitness/fitnes_sum and more chances to be randomly selected
 */
float geneticAlgo::calculateFitnessSum(std::vector<Meal> &mealVector)
{
    float fitness_sum = 0.0;
    float fitness_stage1 =0.0;
    for ( Meal &mealIterator :mealVector)
    {
       fitness_stage1 = mealIterator.calorieInMeal*mealIterator.calorieInMeal;
        if (fitness_stage1 < 1)
        {
            fitness_stage1 = 1;
        }
        mealIterator.fitness_value =  1/fitness_stage1;
        
        fitness_sum = fitness_sum +mealIterator.fitness_value ;
    }
    return fitness_sum;
}
/**
 * @brief Select First N member from meal population
 * 
 * @param mealVector Meal Population to select from
 * @param NumberToSelect 
 */
void geneticAlgo::selectFirstNMeal(std::vector<Meal> &mealVector, const int &NumberToSelect)
{
    mealVector.resize(NumberToSelect);
}

/**
 * @brief Search the parent based on population's commulutive probability
 * 
 * @param mealPopulation Parent population to search from 
 * @param probability randomly generated float value
 * @return int Location of parents on foodData List
 */
int geneticAlgo::binaryProbSearch(const std::vector<Meal> &mealPopulation,const float &probability)
{
    int leftPointer =0;
    int rightPointer = mealPopulation.size()-1;
    int midTerm =0;
    while (leftPointer <= rightPointer)
    {
        midTerm = leftPointer + ((rightPointer - leftPointer)/2);

        //Check if the prob is at mid
        if ( (mealPopulation[midTerm].fitness_cpd[0] < probability) &&  (mealPopulation[midTerm].fitness_cpd[1] > probability ) )
        {
            return midTerm;
        }
        //If probability is less than cpd[0] of mid term than ignore the right half part
        if (probability < mealPopulation[midTerm].fitness_cpd[0])
        {
            rightPointer = midTerm - 1;
        }

        else //mealPopulation[midTerm].fitness_cpd[1] > probability
        {
            leftPointer = midTerm + 1;
        }

    }
    return 0;
}

/**
 * @brief This function calculates fitness value for each meal from meal population
 * 
 * @param mealVector Population of meal whose fitness function has to be calculated
 */
void geneticAlgo::calFitnessProbForEachMeal(std::vector<Meal> &mealVector, const float &fitness_sum)
{
    float lastProb = 0.0;
    for (Meal &mealIterator :mealVector )
    {
        mealIterator.fitness_cpd[0] = lastProb;
        lastProb = lastProb + (mealIterator.fitness_value / fitness_sum);
        mealIterator.fitness_cpd[1] = lastProb;
    }
    _mFitnessSum = 0.0;
}

/**
 * @brief Sort in descending order based on meal fitness_value. The higher the fitness_value, the more required qualities in meal
 * 
 * @param mealVector Meal population to be sorted
 */
void geneticAlgo::fitnessBasedMealStructSort(std::vector<Meal> &mealVector)
{
    std::sort(mealVector.begin(),mealVector.end(),[](const Meal &meal1, const Meal &meal2)
    {
        return meal1.fitness_value > meal2.fitness_value;
    });
}
/**
 * @brief This function generate a random float between 0.0 to 2.0 in steps of 0.1
 * 
 * @return float servingSize for meal
 */
float geneticAlgo::generateRandomServingSize()
{
    float servingSize=0.0;
    auto randomNumberFunc = [](int low, int high)
    {
        auto randomFunc = [distribution_ = std::uniform_int_distribution<int>(low,high),
                            random_engine_ = std::mt19937{std::random_device{}()}]() mutable
                            {
                                return distribution_(random_engine_);
                            };
                            return randomFunc;
    };
    std::vector<int> randomNumber = {-1};
    std::generate(begin(randomNumber),end(randomNumber),randomNumberFunc(0,20));
    servingSize = ((float)randomNumber[0])*_mServingSizeParts; //To get float between 0.0 to 2.0
    //TODO make serving size step as const variable 0.1
    return servingSize;
}
/**
 * @brief This function generate random 5th digit float number
 * 
 * @return float 5th decimal place number to be used as probability
 */
float geneticAlgo::generate5thDigitProb()
{
    float probability =0.0f;
    auto randomNumberFunc = [](int low, int high)
    {
        auto randomFunc = [distribution_ = std::uniform_int_distribution<int>(low,high),
                            random_engine_ = std::mt19937{std::random_device{}()}]() mutable
                            {
                                return distribution_(random_engine_);
                            };
                            return randomFunc;
    };
    std::vector<int> probVec = {0};
    std::generate(begin(probVec),end(probVec),randomNumberFunc(0,100000));
    probability = ((float)probVec[0])*0.00001;
    return probability;
}
/**
 * @brief This generate vector of random integers
 * 
 * @param lastIndex The highest number to select random integer from
 * @param sizeOfVector Number of integers required in vector
 * @return std::vector<int> Vector of random integer of size 'sizeOFVector'
 */
std::vector<int> geneticAlgo::generateRandomIntVector(const int &lastIndex, const int &sizeOfVector )
{
    std::vector<int> randomVec;
    randomVec.reserve(sizeOfVector);
    std::function<std::function<int()>(int,int)>randomNumBetweem = [] (int low,int high)
    {
        std::function<int()> randomFunc = [distribution_ = std::uniform_int_distribution<int>(low,high),
                            random_engine_ = std::mt19937{std::random_device{}()}]() mutable
                            {
                                return distribution_(random_engine_);
                            };
                            return randomFunc;
    };
    std::generate_n(std::back_inserter(randomVec),sizeOfVector,randomNumBetweem(0,lastIndex));
    return randomVec;
}

/**
 * @brief This function generate a random integer   
 * 
 * @param lastIndex The highest number to select random integer from
 * @return int random Integer
 */
int geneticAlgo::generateRandomInt(const int &lastIndex)
{
    auto randomNumberFunc = [](int low, int high)
    {
        auto randomFunc = [distribution_ = std::uniform_int_distribution<int>(low,high),
                            random_engine_ = std::mt19937{std::random_device{}()}]() mutable
                            {
                                return distribution_(random_engine_);
                            };
                            return randomFunc;
    };
    std::vector<int> randomNumber = {-1};
    std::generate(begin(randomNumber),end(randomNumber),randomNumberFunc(0,lastIndex));
    return randomNumber[0];
}
