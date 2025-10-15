"""
Folder management utilities for pipeline execution.

Handles automatic creation and organization of all pipeline-related directories
following the established conventions from the Jupyter notebooks.
"""

import os
from typing import Dict

try:
    from .config_manager import ConfigManager
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from config_manager import ConfigManager


class FolderManager:
    """
    Manages folder structure for pipeline execution.
    
    Creates and maintains the standard directory hierarchy used by all tools,
    following conventions established in the Jupyter notebooks.
    """
    
    def __init__(self, pipeline_name: str, job_folder: str, debug: bool=False, shared: bool=False):
        """
        Initialize folder manager for a pipeline.

        Args:
            pipeline_name: Name of the pipeline (used for output folders)
            job_folder: Name of the job (used for output folders)
            debug: if True, it will not attempt to generate folder. Useful to check locally the expected input output of models
            shared: if True, uses 'public' as the user name for shared/public pipelines
        """
        self._folders = {}

        # Load configuration
        config_manager = ConfigManager()
        folder_config = config_manager.get_folder_config()

        """Setup base folder structure following notebook conventions."""
        # Get user info from current directory structure
        notebooks_folder = os.getcwd() if not debug else "biopipelines"
        user_name = os.path.basename(os.path.dirname(notebooks_folder)) if not debug else "USER"
        self.user_name = user_name
        user_code = user_name if not shared else "public"

        # Get base paths from config
        group_folder = folder_config.get('group', '/shares/locbp.chem.uzh')

        # Base system folders
        # For shared pipelines, use a shared PDBs folder accessible to all users
        if shared:
            pdbs_folder = os.path.join(group_folder, user_code, "PDBs")
        else:
            pdbs_folder = os.path.join(notebooks_folder, "PDBs")

        self._folders.update({
            "notebooks": notebooks_folder,
            "home": f"/home/{user_name}",
            "data": f"/data/{user_name}",
            "scratch": f"/scratch/{user_name}",
            "group": group_folder,
            "user": f"{group_folder}/{user_code}",
            "PDBs": pdbs_folder,
            "HelpScripts": os.path.join(notebooks_folder, "HelpScripts"),
            "MMseqs2": os.path.join(notebooks_folder, "MMseqs2")})

        # Shared folders from config
        shared_config = folder_config.get('shared', {})
        models_path = shared_config.get('models', 'models')
        containers_path = shared_config.get('containers', 'containers')
        boltz_cache_path = shared_config.get('BoltzCache', 'models/Boltz')

        self._folders.update({
            "models": os.path.join(self._folders["group"], models_path),
            "containers": os.path.join(self._folders["group"], containers_path)})

        # Tool-specific data folders from config
        tool_data_config = folder_config.get('tool_data', {})
        dynamicbind_base = os.path.join(self._folders["data"], tool_data_config.get('DynamicBind', 'DynamicBind'))
        self._folders.update({
            "RFdiffusion": os.path.join(self._folders["data"], tool_data_config.get('RFdiffusion', 'RFdiffusion')),
            "RFdiffusionAllAtom": os.path.join(self._folders["data"], tool_data_config.get('RFdiffusionAllAtom', 'rf_diffusion_all_atom')),
            "ProteinMPNN": os.path.join(self._folders["data"], tool_data_config.get('ProteinMPNN', 'ProteinMPNN')),
            "AlphaFold": os.path.join(self._folders["data"], tool_data_config.get('AlphaFold', 'localcolabfold')),
            "OmegaFold": os.path.join(self._folders["data"], tool_data_config.get('OmegaFold', 'OmegaFold')),
            "Boltz": os.path.join(self._folders["data"], tool_data_config.get('Boltz', 'boltz')),
            "BoltzCache": os.path.join(self._folders["group"], boltz_cache_path),
            "DynamicBind": dynamicbind_base,
            "DynamicBindWeights": os.path.join(dynamicbind_base, "workdir"),
            "BioPipelines": os.path.join(self._folders["user"], "BioPipelines")})
        self._folders.update({
            "pipeline": os.path.join(self._folders["BioPipelines"],pipeline_name) if not debug else "Debug"})

        for folder_key in ["user", "BioPipelines", "PDBs"]:
            if not debug: os.makedirs(self._folders[folder_key], exist_ok=True) #create directories
        os.makedirs(self._folders["pipeline"], exist_ok=True)

        self._folders["output"] = FolderManager.unique_name(self._folders["pipeline"],job_folder,full_path=True)
        os.makedirs(self._folders["output"], exist_ok=True)
        
        # Extract job ID from the output folder name
        output_basename = os.path.basename(self._folders["output"])
        if "_" in output_basename:
            self.job_id = output_basename.split("_")[-1]
        else:
            raise ValueError(f"Could not extract job ID from output folder name: {output_basename}")
        self._folders["runtime"] = os.path.join(self._folders["output"], "RunTime")
        os.makedirs(self._folders["runtime"], exist_ok=True)

        self._folders["logs"] = os.path.join(self._folders["output"], "Logs")
        os.makedirs(self._folders["logs"], exist_ok=True)
        
    def get_folders(self) -> Dict[str, str]:
        """
        Get dictionary of all folder paths.
        
        Returns:
            Dictionary mapping folder names to absolute paths
        """
        return self._folders.copy()
    
    def unique_name(directory: str, root: str, ext: str = "", w: int = 3, full_path: bool = False) -> str:
        """
        Generate unique filename/foldername to avoid conflicts.
        
        Based on the unique_name function from notebooks.
        
        Args:
            directory: Directory to check for existing files
            root: Base name for file/folder
            ext: File extension (optional)
            w: Width of numeric suffix (default 3)
            
        Returns:
            Unique name that doesn't exist in directory
        """
        i = 1
        while True:
            u_name = root + "_" + "{:0>{width}}".format(i, width=w)
            if not os.path.exists(os.path.join(directory, u_name + ext)):
                if full_path:   return os.path.join(directory, u_name + ext)
                else:           return u_name + ext
            i += 1
    
    def __getitem__(self, key: str) -> str:
        """Allow dictionary-like access to folders."""
        return self._folders[key]
    
    def __contains__(self, key: str) -> bool:
        """Check if folder key exists."""
        return key in self._folders