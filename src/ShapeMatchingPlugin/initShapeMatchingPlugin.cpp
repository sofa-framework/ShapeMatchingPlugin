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

namespace sofa::component
{

	extern "C" {
		SOFA_SHAPEMATCHINGPLUGIN_API void initExternalModule();
		SOFA_SHAPEMATCHINGPLUGIN_API const char* getModuleName();
		SOFA_SHAPEMATCHINGPLUGIN_API const char* getModuleVersion();
		SOFA_SHAPEMATCHINGPLUGIN_API const char* getModuleLicense();
		SOFA_SHAPEMATCHINGPLUGIN_API const char* getModuleDescription();
		SOFA_SHAPEMATCHINGPLUGIN_API const char* getModuleComponentList();
	}

	void initShapeMatchingPlugin()
	{
		initExternalModule();
	}

	void initExternalModule()
	{
		static bool first = true;
		if (first)
		{
			first = false;
		}
	}

	const char* getModuleName()
	{
		return shapematching::MODULE_NAME;
	}

	const char* getModuleVersion()
	{
		return shapematching::MODULE_VERSION;
	}

	const char* getModuleLicense()
	{
		return "Private";
	}

	const char* getModuleDescription()
	{
		return "Plugin with ShapeMatchingPlugin";
	}

	const char* getModuleComponentList()
	{
		/// string containing the names of the classes provided by the plugin
		static std::string classes = sofa::core::ObjectFactory::getInstance()->listClassesFromTarget(sofa_tostring(SOFA_TARGET));
		return classes.c_str();
	}

} // namespace sofa::component
