/*
 * JSONTest.cpp
 *
 *  Created on: Nov 6, 2015
 *      Author: pfandedd
 */

#include "../src/sgpp/base/tools/JSON.hpp"

int main(int argc, char **argv) {

  SGPP::base::JSON configuration("test.json");

//  SGPP::base::JSONNode &node = configuration["Inhaber"];
//  SGPP::base::JSONNode &subNode = node["Name"];

  std::cout << "Herausgeber -> " << configuration["Herausgeber"].getValue() << std::endl;

  std::cout << "Inhaber/Name -> " << configuration["Inhaber"]["Name"].getValue() << std::endl;

  SGPP::base::JSONNode &listNode = configuration["Inhaber"]["Hobbys"];

  std::cout << "Inhaber/Hobbies -> ";
  bool first = true;
  for (size_t i = 0; i < listNode.size(); i++) {
    SGPP::base::JSONNode &valueNode = listNode.getItem(i);
    if (first) {
      first = false;
    } else {
      std::cout << ", ";
    }
    std::cout << valueNode.getValue();
  }

  std::cout << std::endl;

  configuration.serialize("write.json");

  return 0;
}
