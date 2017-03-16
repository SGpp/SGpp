// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>

#include <sgpp/base/tools/json/TextNode.hpp>

int main(int argc, char** argv) {
  try {
    json::JSON configuration;

    configuration.addDictAttr("parent").addTextAttr("t1", "v1").addTextAttr("t2",
        "v2").addListAttr("list1").addTextValue(
          "tv1").addIdValue(96.0).addTextValue("tv2");
    configuration.addTextAttr("textAttr1", "text1").addIDAttr("numVal1", 36.0);

    std::cout << "value: " << configuration["parent"]["list1"][1].get() <<
              std::endl;

    configuration["parent"]["list1"][1].setDouble(7);

    std::cout << "value: " << configuration["parent"]["list1"][1].get() <<
              std::endl;


    std::cout << "value: " << configuration["parent"]["list1"][1].getDouble() <<
              std::endl;

    configuration["numVal1"].erase();

    std::unique_ptr<json::Node> parentNode = configuration["parent"].erase();

    configuration.addDictAttr("parentparent").addAttribute("parent",
        std::move(parentNode));

    configuration.serialize("write.json");

    json::JSON reread("write.json");

    reread.serialize("write2.json");
  } catch (json::json_exception& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
