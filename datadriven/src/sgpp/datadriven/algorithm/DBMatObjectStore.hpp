// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflinePermutable.hpp>
#include <sgpp/datadriven/configuration/GeometryConfiguration.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

class DBMatObjectStore {
 public:
  /**
   * @brief Default constructor.
   *
   */
  DBMatObjectStore();

  /**
   * @brief Constructor with path to database file.
   *
   * @param fileName
   */
  explicit DBMatObjectStore(const std::string& fileName);

  /**
   * @brief Stores a given offline object together with its configuration in the object store.
   *
   * @param gridConfig Grid configuration
   * @param geometryConfig Geometry configuration for geometry aware sparse grids
   * @param adaptConfig Adaptivity configuration
   * @param regularizationConfig Regularization configuration
   * @param densityEstimationConfig Density estimation configuration
   * @param object The object to be stored
   */
  void putObject(const sgpp::base::GeneralGridConfiguration& gridConfig,
                 const sgpp::datadriven::GeometryConfiguration& geometryConfig,
                 const sgpp::base::AdaptivityConfiguration& adaptConfig,
                 const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                 const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
                 const DBMatOffline* object);
  /**
   * @brief Returns a suitable base object for the permutation and blow-up approach.
   * The grid configuration of the stored base object is returned in baseGridConfig.
   * If no suitable base object exists, a nullptr is returned.
   *
   * @param gridConfig Grid configuration of the desired offline object
   * @param geometryConfig Geometry configuration for geometry aware sparse grids
   * @param adaptConfig Adaptivity configuration
   * @param regularizationConfig Regularization configuration
   * @param densityEstimationConfig Density estimation configuration
   * @param baseGridConfig Reference to a grid configuration. Gets overridden by the grid
   * configuration of the returned base object
   * @return const DBMatOfflinePermutable*
   */
  const DBMatOfflinePermutable* getBaseObject(
      const sgpp::base::GeneralGridConfiguration& gridConfig,
      const sgpp::datadriven::GeometryConfiguration& geometryConfig,
      const sgpp::base::AdaptivityConfiguration& adaptConfig,
      const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
      sgpp::base::GeneralGridConfiguration& baseGridConfig);

  /**
   * @brief Returns an identical offline object to the specified configuration.
   * If no such object exits, a nullptr is returned
   *
   * @param gridConfig Grid configuration
   * @param geometryConfig Geometry configuration for geometry aware sparse grids
   * @param adaptConfig Adaptivity configuration
   * @param regularizationConfig Regularization configuration
   * @param densityEstimationConfig Density estimation configuration
   * @return const DBMatOffline*
   */
  const DBMatOffline* getObject(
      const sgpp::base::GeneralGridConfiguration& gridConfig,
      const sgpp::datadriven::GeometryConfiguration& geometryConfig,
      const sgpp::base::AdaptivityConfiguration& adaptConfig,
      const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig);

 protected:
  /**
   * @brief Datastructure to store offline objects together with their configuration.
   * This class is not intended to be used outside of DBMatObjectStore.
   */
  class ObjectContainer {
   public:
    /**
     * @brief Public constructor. Gets initialized with the offline object and the correspoding
     * configuration. Note that ownership of the given offline object get transfered to the
     * container. I.e. if the container is deleted, its offline object is deleted as well.
     *
     * @param gridConfig Grid configuration
     * @param geometryConfig Geometry configuration for geometry aware sparse grids
     * @param adaptConfig Adaptivity configuration
     * @param regularizationConfig Regularization configuration
     * @param densityEstimationConfig Density estimation configuration
     * @param offlineObject Unique pointer to an offline object. Ownership of this object gets
     * transfered to the container
     */
    explicit ObjectContainer(
        const sgpp::base::GeneralGridConfiguration& gridConfig,
        const sgpp::datadriven::GeometryConfiguration& geometryConfig,
        const sgpp::base::AdaptivityConfiguration& adaptConfig,
        const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
        const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
        std::unique_ptr<const DBMatOffline> offlineObject);

    /**
     * @brief Returns a read-only reference to the containers offline object.
     *
     * @return const DBMatOffline&
     */
    const DBMatOffline& getOfflineObject() const;

    /**
     * @brief Returns a read-only reference to the containers grid configuration.
     *
     * @return const sgpp::base::GeneralGridConfiguration&
     */
    const sgpp::base::GeneralGridConfiguration& getGridConfig() const;

    /**
     * @brief Checks wheter the configuration of a container matches a given configuration.
     * If searcBase = true, it is checked wheter the offline object is a suitable base object for
     * the permutation and blow-up approach.
     *
     * @param gridConfig Grid configuration
     * @param geometryConfig Geometry configuration for geometry aware sparse grids
     * @param adaptConfig Adaptivity configuration
     * @param regularizationConfig Regularization configuration
     * @param densityEstimationConfig Density estimation configuration
     * @param searchBase Flag to specify whether an identical offline object or a suitable base
     * object is to be searched
     * @return true
     * @return false
     */
    bool configMatches(
        const sgpp::base::GeneralGridConfiguration& gridConfig,
        const sgpp::datadriven::GeometryConfiguration& geometryConfig,
        const sgpp::base::AdaptivityConfiguration& adaptConfig,
        const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
        const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
        bool searchBase);

   private:
    // Member variables, i.e. the offline object and its configuration.
    sgpp::base::GeneralGridConfiguration gridConfig;
    sgpp::datadriven::GeometryConfiguration geometryConfig;
    sgpp::base::AdaptivityConfiguration adaptConfig;
    RegularizationConfiguration regularizationConfig;
    DensityEstimationConfiguration densityEstimationConfig;
    std::unique_ptr<const DBMatOffline> offlineObject;
  };
  // Vector with stored object containers.
  std::vector<ObjectContainer> objects;
  // Optional path to a database file
  std::string dbFilePath;
  // True if database file is given
  bool hasDatabase;

  /**
   * @brief Returns the index to a suitable offline object.
   * If searchBase = true, a suitable base object for the permutation and blow-up approach is
   * searched for. If no suitable object exists, SIZE_MAX is returned.
   *
   * @param gridConfig Grid configuration
   * @param geometryConfig Geometry configuration for geometry aware sparse grids
   * @param adaptConfig Adaptivity configuration
   * @param regularizationConfig Regularization configuration
   * @param densityEstimationConfig Density estimation configuration
   * @param searchBase Flag to specify whether an identical offline object or a suitable base object
   * is to be searched
   * @return int
   */
  size_t getObjectContainerIndex(
      const sgpp::base::GeneralGridConfiguration& gridConfig,
      const sgpp::datadriven::GeometryConfiguration& geometryConfig,
      const sgpp::base::AdaptivityConfiguration& adaptConfig,
      const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
      bool searchBase = false);

  /**
   * @brief Returns the object container for the given index.
   *
   * @param index Index of the object container
   * @return const ObjectContainer&
   */
  const ObjectContainer& getObjectContainer(size_t index) const;
};

}  // namespace datadriven
}  // namespace sgpp
