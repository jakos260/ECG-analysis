# ECG-analysis - PhD project on modeling heart repolarization

This repository contains the codebase and analytical tools for my PhD project focused on the modeling of heart repolarization. It includes data processing algorithms, visualization helpers, and the core analysis scripts written in MATLAB and Python.

## 📂 Repository Structure

The repository is organized into the following main directories:

```text
.
├── data/                  # ⚠️ MUST BE POPULATED MANUALLY (see Setup below)
├── Matlab/                # MATLAB scripts and functions
│   ├── ecg analysis/      # Main PhD project directory (core repolarization models)
│   └── ...                # Data processing helpers and additional functions
└── Python/                # Python utilities
    ├── qtripy             # MATLAB to QtiPlot communication scripts
    └── ...                # Data visualization helpers
```

## ⚙️ Environment Setup

### MATLAB Configuration

To run the MATLAB scripts correctly without tracking local paths in version control, you need to set up an untracked local configuration file.

1. **Create a local paths file:**
   Navigate to the `Matlab/ecg_analysis/` directory and create a new file named `local_paths.m` (this file is intentionally ignored by Git).

2. **Define your local data directory:**
   Add the following line to `local_paths.m`, replacing the path with the absolute path to your local `data` directory:
   ```matlab
   % local_paths.m
   dir_name = 'C:\path\to\your\Data';
   ```
   
The `startup.m` script is already configured to automatically detect and run this file, injecting the `dir_name` variable into your workspace.

### Python Configuration

To use the Python utilities, it is recommended to set up an isolated virtual environment. Follow these steps:

1. **Navigate to the Python directory:**
   ```bash
   cd Python
   ```

2. **Create a virtual environment:**
   ```bash
   python -m venv ecg_project_venv
   ```

3. **Activate the environment:**
   * **Windows:**
     ```cmd
     ecg_project_venv\Scripts\activate
     ```
   * **macOS / Linux:**
     ```bash
     source ecg_project_venv/bin/activate
     ```

4. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

5. **Environment Variables:**

The project relies on specific absolute paths to locate data and external software. You need to create a `.env` file in your working directory and define the following variables:

```env
# Define absolute paths for your local setup
ENV_DATA_PATH="...\Data"
ENV_QTRIPLOT_PATH="...\qtriplot.exe"
```
