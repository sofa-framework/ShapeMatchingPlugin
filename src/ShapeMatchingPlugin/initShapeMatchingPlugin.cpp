/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <ShapeMatchingPlugin/initShapeMatchingPlugin.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/system/PluginManager.h>

namespace shapematchingplugin
{

	extern void registerShapeMatchingForceField(sofa::core::ObjectFactory* factory);
	extern void registerShapeMatchingRotationFinder(sofa::core::ObjectFactory* factory);

	extern "C" {
		SOFA_SHAPEMATCHINGPLUGIN_API void initExternalModule();
		SOFA_SHAPEMATCHINGPLUGIN_API const char* getModuleName();
		SOFA_SHAPEMATCHINGPLUGIN_API const char* getModuleVersion();
		SOFA_SHAPEMATCHINGPLUGIN_API const char* getModuleLicense();
		SOFA_SHAPEMATCHINGPLUGIN_API const char* getModuleDescription();
		SOFA_SHAPEMATCHINGPLUGIN_API void registerObjects(sofa::core::ObjectFactory* factory);
	}

	void initShapeMatchingPlugin()
	{
		initExternalModule();
	}

	void initExternalModule()
	{
		// make sure that this plugin is registered into the PluginManager
        sofa::helper::system::PluginManager::getInstance().registerPlugin(MODULE_NAME);

		static bool first = true;
		if (first)
		{
			first = false;
		}
	}

	const char* getModuleName()
	{
		return MODULE_NAME;
	}

	const char* getModuleVersion()
	{
		return MODULE_VERSION;
	}

	const char* getModuleLicense()
	{
		return "Private";
	}

	const char* getModuleDescription()
	{
		return "Plugin with ShapeMatchingPlugin";
	}

	void registerObjects(sofa::core::ObjectFactory* factory)
	{
		void registerShapeMatchingForceField(factory);
		void registerShapeMatchingRotationFinder(factory);
	}

} // namespace shapematchingplugin
