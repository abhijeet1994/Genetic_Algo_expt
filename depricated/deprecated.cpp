template <typename A, typename B>
    void zip (const std::vector<A> &a, const std::vector<B> &b, std::vector<std::pair<A,B>> &zipped);
    template <typename A, typename B>
    void unzip (const std::vector<std::pair<A,B>> &zipped, std::vector<A> &a, std::vector<B> &b);
template <typename A, typename B>
void geneticAlgo::zip(const std::vector<A> &a, const std::vector<B> &b, std::vector<std::pair<A,B>> &zipped)
{
    for (size_t i=0; i < a.size(); ++i)
    {
        zipped.push_back(std::make_pair(a[i],b[i]));
    }
}
template <typename A, typename B>
void geneticAlgo::unzip(const std::vector<std::pair<A,B>> &zipped, std::vector<A> &a, std::vector<B> &b)
{
    for (size_t i=0; i < a.size(); ++i)
    {
        a[i] = zipped[i].first;
        b[i] = zipped[i].second;
    }
}
void calorieBasedSort(std::vector<std::vector<int>> &mealItems,std::vector<float> &calorieInMeal);
void geneticAlgo::calorieBasedSort(std::vector<std::vector<int>> &mealItems,std::vector<float> &calorieInMeal)
{
    std::cout <<"before state of Main vector " << std::endl;
      for (size_t i=0; i< calorieInMeal.size(); ++i)
    {
        std::cout << "Calorie"  << calorieInMeal[i] <<" ' " <<std::endl;
        // std::cout <<  std::endl;
        std::cout << "MEalItem " <<mealItems[i][0] <<  mealItems[i][1] <<  mealItems[i][2] << std::endl;
        if (i==3) 
        {
            break;
        }
    }
    std::cout << std::endl;
    std::vector<std::pair<std::vector<int>,float>> zipped; 
    zip(mealItems,calorieInMeal,zipped);
    std::sort(std::begin(zipped), std::end(zipped),
    [&](const auto &a, const auto &b)
    {
        return a.second < b.second;
    }
    );
    unzip(zipped,mealItems,calorieInMeal);
    std::cout <<"After State of CalorieVec " <<std::endl;
    for (size_t i=0; i< calorieInMeal.size(); ++i)
    {
        std::cout <<" calorie, "<< calorieInMeal[i] << "," ;
        std::cout <<  std::endl;
        std::cout <<" mealItem, " << mealItems[i][0] <<  mealItems[i][1] <<  mealItems[i][2] << std::endl;
        if (i==3) 
        {
            break;
        }
    }
    std::cout <<std::endl;
}
    void subtractIdealCalorieFromMeal(const float &idealCalorie);

void geneticAlgo::subtractIdealCalorieFromMeal(const float &idealCalorie)
{
    for (float iterator : _mCalorieInMeal)
    {
        iterator = iterator - idealCalorie;
    }
}

void geneticAlgo::calculateCalorieInMeal()
{
    float calorie =0.0;
    float protein_cont = 0.0;
    float carbs_cont = 0.0;
    float fats_cont = 0.0;
    std::vector<float> macro;
    macro.reserve(3);
    for (size_t counter=0; counter<_mPopulationSize;counter++)
    {
        for (int j=0; j < _mNumCompoInRule; j++)
        {
            calorie = calorie + ((_mfoodList[_mCompPositionInRule[j]][_mMealItems[counter][j]].kiloCalories)*_mServingSize[counter][j]);
            protein_cont = protein_cont + ((_mfoodList[_mCompPositionInRule[j]][_mMealItems[counter][j]].proteinContent)*_mServingSize[counter][j]);
            macro.push_back(protein_cont);
            carbs_cont = carbs_cont + ((_mfoodList[_mCompPositionInRule[j]][_mMealItems[counter][j]].carbsContent)*_mServingSize[counter][j]);
            macro.push_back(carbs_cont);
            fats_cont = fats_cont + ((_mfoodList[_mCompPositionInRule[j]][_mMealItems[counter][j]].FatContent)*_mServingSize[counter][j]);
            macro.push_back(fats_cont);
        }
        _mCalorieInMeal.push_back(calorie);
        _mMacrosInMeal.push_back(macro);
        macro.clear();
    }
}
    void calculateCalorieInMeal();

    void createPopulationOfMeal();

}
/**
 * @brief This creates meal population as a vector 
 * 
 */
void geneticAlgo::createPopulationOfMeal()
{
    std::vector<int> individualMeal;
    std::vector<float> individualServeSize;
    _mMealItems.reserve(_mPopulationSize);
    _mServingSize.reserve(_mPopulationSize);
    
    for (size_t i =0; i < _mPopulationSize ; i++)
    {
        individualMeal.clear();
        individualServeSize.clear();
        individualMeal.reserve(_mNumCompoInRule);
        individualServeSize.reserve(_mNumCompoInRule);
        for (int j=0; j< _mNumCompoInRule; j++)
        {
            individualMeal.push_back(generateRandomInt(_mfoodList[_mCompPositionInRule[j]].size()));
            individualServeSize.push_back(generateRandomInt(_mfoodList[_mCompPositionInRule[j]].size()));
        }
        _mMealItems.push_back(individualMeal);
        _mServingSize.push_back(individualServeSize);
    }
}

