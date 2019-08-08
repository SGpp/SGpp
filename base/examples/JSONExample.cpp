// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_json_cpp Using JSON
 *
 * This example demonstrates how to use the basic functionality of SG++ JSON API.
 * The goal of this API is to facilitate generation and parsing of complex configuration files.
 * Also it permits easy but limited object serialization functionality.
 */

#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/TextNode.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>

int main(int argc, char **argv) {
  try {
    /**
     * We first create an empty JSON object
     */
    json::JSON configuration;
    /**
     * Then we introduce the following JSON structure:
     * @verbatim
     * {
     * 	parent:{
     * 		"t1":"v1";
     * 		"t2":"v2";
     * 		"list1":["tv1",96.0,"tv2"]
     * 			};
     * 	"textAttr1":"text1";
     * 	"numVal1":36.0
     * }
     * @endverbatim
     *
     * As you can see in the code below,
     * methods to append new entries support chaining for convenience.
     */
    configuration.addDictAttr("parent")
        .addTextAttr("t1", "v1")
        .addTextAttr("t2", "v2")
        .addListAttr("list1")
        .addTextValue("tv1")
        .addIdValue(96.0)
        .addTextValue("tv2");
    configuration.addTextAttr("textAttr1", "text1").addIDAttr("numVal1", 36.0);

    /**
     * You can conveniently get and set values inside the JSON object using Key/Value notation.
     *
     * This call gets the string representation of the second entry (<tt>96.0</tt>) in
     * <tt>list1</tt> of
     * dictionary <tt>parent</tt>.
     */
    std::cout << "value: " << configuration["parent"]["list1"][1].get() << std::endl;

    /**
     * This call sets the second entry (<tt>96.0</tt>) in <tt>list1</tt> of
     * dictionary <tt>parent</tt> to the double value <tt>7</tt> and then prints its string
     * representation.
     */
    configuration["parent"]["list1"][1].setDouble(7);
    std::cout << "value: " << configuration["parent"]["list1"][1].get() << std::endl;

    /**
     * You can also try to obtain a typed representation of an entry. Conversion is performed
     * automatically and may result in  an exception if it can't be converted.
     */
    std::cout << "value: " << configuration["parent"]["list1"][1].getDouble() << std::endl;

    /**
     * Next we demonstrate erasing values.
     *
     * Here we erase <tt>numVal1</tt> from the top level of the <tt>configuration</tt> object.
     */
    configuration["numVal1"].erase();

    /**
     * Now we erase the dictionary named <tt>parent</tt> from the the <tt>configuration</tt> object.
     * Note
     * that the object that is erased is
     * returned to us by the erase operation. All Objects inside a JSON object are represented
     * internally as an object of type <tt>Node</tt>"
     */
    std::unique_ptr<json::Node> parentNode = configuration["parent"].erase();

    /**
     *
     */
    configuration.addDictAttr("parentparent").addAttribute("parent", std::move(parentNode));

    /**
     * Serialization (i.e. writing of a standard compatible json file) is also possible.
     */
    configuration.serialize("write.json");

    /**
     * We can also read a JSON file by passing it's path to the constructor of a JSON object
     */
    json::JSON reread("write.json");
  } catch (json::json_exception &e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
