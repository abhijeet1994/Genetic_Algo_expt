#include "getFoodData.h"

getFoodData::getFoodData(std::string &aFileLocation)
:_mFileLocation(aFileLocation),
_mDlimeter(",")
{
    _mFoodData.foodCategories.clear();
    _mFoodData.food_database.clear();
    _mFoodData.food_database.reserve(foodDataConsts::NumOfCategories);
}
/**
 * @brief This function will convert csv data into internal c++ structure
 * 
 */

void getFoodData::getCSVdata()
{
    std::ifstream file(_mFileLocation);
    std::string line = "";
    std::vector<foodData> food_category;
    food_category.reserve(foodDataConsts::varitiesOfFood);
    std::vector<std::string> vec;
    vec.reserve(foodDataConsts::NumOfColumInFoodData);
    // Iterate through each line and split the content using delimeter
    int num_categories = 0;
    while (getline(file, line))
    {
        vec.clear();   
        boost::algorithm::split(vec, line, boost::is_any_of(_mDlimeter));
        if (_mFoodData.foodCategories.find(vec[2]) == _mFoodData.foodCategories.end()){
            _mFoodData.foodCategories[vec[2]] = num_categories;
            _mFoodData.food_database.push_back(food_category);
            num_categories += 1;
        }
        _mFoodData.food_database[_mFoodData.foodCategories[vec[2]]].push_back(storeDataInVector(vec));
    }
    // Close the File
    file.close();
}
    

/**
 * @brief This function save component of food into its foodData structure, 
 * append foodData to respective food type vector 
 * 
 * @param foodDataVector Vector that hold food compoenet according to its type 
 * @param inputvector food information from csv file in form of vector of strings
 */

foodData getFoodData::storeDataInVector(std::vector<std::string> &inputvector)
{
    foodData data;
    data.serialNumber = ::atoi(inputvector[0].c_str());
    data.ingredientName = inputvector[1];
    data.classification = inputvector[2];
    data.type = inputvector[3];
    data.proteinContent = ::atof(inputvector[4].c_str());
    data.carbsContent = ::atof(inputvector[5].c_str());
    data.FatContent = ::atof(inputvector[6].c_str());
    data.kiloCalories = ::atof(inputvector[7].c_str());
    data.price = ::atof(inputvector[8].c_str());
    data.quantity = ::atof(inputvector[9].c_str());
    return data;
}

/**
 * @brief This function returns mega list which combine all food item
 * 
 * @return std::vector<std::vector<foodData>>  List of all food types
 */
FoodDatabase getFoodData::getFoodDataList()
{
    getCSVdata();
    return _mFoodData;
}