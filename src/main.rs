use eframe::egui;
use rfd::FileDialog;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::process::{Child, Command};
use std::time::{Duration, Instant};
use std::{default, fs};

#[derive(Serialize, Deserialize, Debug, Default)]
struct EnzymeConfig {
    missed_cleavages: i32,
    min_len: i32,
    max_len: i32,
    cleave_at: String,
    restrict: String,
    c_terminal: bool,
}

#[derive(Serialize, Deserialize, Debug, Default)]
struct DatabaseConfig {
    bucket_size: i32,
    enzyme: EnzymeConfig,
    fragment_min_mz: f64,
    fragment_max_mz: f64,
    peptide_min_mass: f64,
    peptide_max_mass: f64,
    ion_kinds: Vec<String>,
    min_ion_index: i32,
    max_variable_mods: i32,
    decoy_tag: String,
    generate_decoys: bool,
    fasta: String,
}

#[derive(Serialize, Deserialize, Debug)]
enum ToleranceConfig {
    #[serde(rename = "da")]
    Da(f64, f64),
    #[serde(rename = "ppm")]
    Ppm(f64, f64),
}

impl Default for ToleranceConfig {
    fn default() -> Self {
        Self::Ppm(-10.0, 10.0)
    }
}

#[derive(Serialize, Deserialize, Debug, Default)]
struct Config {
    database: DatabaseConfig,
    precursor_tol: ToleranceConfig,
    fragment_tol: ToleranceConfig,
    precursor_charge: Vec<i32>,
    isotope_errors: Vec<i32>,
    deisotope: bool,
    chimera: bool,
    wide_window: bool,
    predict_rt: bool,
    min_peaks: i32,
    max_peaks: i32,
    min_matched_peaks: i32,
    max_fragment_charge: i32,
    report_psms: i32,
    mzml_paths: Vec<String>,
}

struct SageLauncher {
    config: Config,
    sage_executable_path: String,
    status_message: String,
    precursor_tolerance_type: ToleranceType,
    fragment_tolerance_type: ToleranceType,
    temp_mzml_path: String,
    process: Option<(Child, Instant)>, // Store both process handle and start time
    elapsed_time: String,
}

#[derive(PartialEq)]
enum ToleranceType {
    PPM,
    DA,
}

impl ToleranceType {
    fn get_default_tolerance(&self) -> ToleranceConfig {
        match self {
            ToleranceType::PPM => ToleranceConfig::Ppm(-10.0, 10.0),
            ToleranceType::DA => ToleranceConfig::Da(-0.5, 0.5),
        }
    }

    fn other(&self) -> ToleranceType {
        match self {
            ToleranceType::PPM => ToleranceType::DA,
            ToleranceType::DA => ToleranceType::PPM,
        }
    }
}

impl Default for SageLauncher {
    fn default() -> Self {
        Self {
            sage_executable_path: String::new(),
            config: Config {
                database: DatabaseConfig {
                    bucket_size: 32768,
                    enzyme: EnzymeConfig {
                        missed_cleavages: 2,
                        min_len: 5,
                        max_len: 50,
                        cleave_at: "KR".to_string(),
                        restrict: "P".to_string(),
                        c_terminal: true,
                    },
                    fragment_min_mz: 200.0,
                    fragment_max_mz: 2000.0,
                    peptide_min_mass: 500.0,
                    peptide_max_mass: 5000.0,
                    ion_kinds: vec!["b".to_string(), "y".to_string()],
                    min_ion_index: 2,
                    max_variable_mods: 2,
                    decoy_tag: "rev_".to_string(),
                    generate_decoys: false,
                    fasta: String::new(),
                },
                precursor_tol: ToleranceConfig::default(),
                fragment_tol: ToleranceConfig::default(),
                precursor_charge: vec![2, 4],
                isotope_errors: vec![-1, 3],
                deisotope: false,
                chimera: false,
                wide_window: false,
                predict_rt: false,
                min_peaks: 15,
                max_peaks: 150,
                min_matched_peaks: 6,
                max_fragment_charge: 1,
                report_psms: 1,
                mzml_paths: Vec::new(),
            },
            status_message: String::new(),
            precursor_tolerance_type: ToleranceType::PPM,
            fragment_tolerance_type: ToleranceType::PPM,
            temp_mzml_path: String::new(),
            process: None,
            elapsed_time: String::new(),
        }
    }
}

impl eframe::App for SageLauncher {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Update process status and elapsed time
        self.update_process_status();
        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.heading("Sage Launcher");

                // Process Status Section
                if let Some((process, start_time)) = &self.process {
                    ui.add_space(10.0);
                    ui.horizontal(|ui| {
                        let status_color = egui::Color32::GREEN;
                        ui.colored_label(status_color, "●"); // Status dot
                        ui.label("Process Running");
                        ui.label(format!("({})", self.elapsed_time));

                        if ui.button("Stop Process").clicked() {
                            self.stop_process();
                        }
                    });
                    ui.add_space(10.0);
                }
                ui.add_space(20.0);

                // File Selection Section
                ui.collapsing("File Selection", |ui| {
                    // Executable file picker
                    ui.horizontal(|ui| {
                        ui.label("Sage Executable:");
                        ui.text_edit_singleline(&mut self.sage_executable_path);
                        if ui.button("Browse").clicked() {
                            if let Some(path) = FileDialog::new().pick_file() {
                                self.sage_executable_path = path.display().to_string();
                            }
                        }
                    });

                    // FASTA file picker
                    ui.horizontal(|ui| {
                        ui.label("FASTA File:");
                        ui.text_edit_singleline(&mut self.config.database.fasta);
                        if ui.button("Browse").clicked() {
                            if let Some(path) = FileDialog::new()
                                .add_filter("FASTA", &["fasta"])
                                .pick_file()
                            {
                                self.config.database.fasta = path.display().to_string();
                            }
                        }
                    });

                    // mzML file picker
                    ui.horizontal(|ui| {
                        ui.label("mzML File:");
                        ui.text_edit_singleline(&mut self.temp_mzml_path);
                        if ui.button("Browse").clicked() {
                            if let Some(path) = FileDialog::new()
                                .add_filter("mzML", &["mzML", "gz", "mzml"])
                                .pick_file()
                            {
                                self.temp_mzml_path = path.display().to_string();
                                self.config.mzml_paths = vec![self.temp_mzml_path.clone()];
                            }
                        }
                    });
                });

                // Database Configuration Section
                ui.collapsing("Database Configuration", |ui| {
                    ui.add(
                        egui::Slider::new(&mut self.config.database.bucket_size, 8192..=65536)
                            .text("Bucket Size"),
                    );

                    // Enzyme Configuration
                    ui.group(|ui| {
                        ui.heading("Enzyme Settings");
                        ui.add(
                            egui::Slider::new(
                                &mut self.config.database.enzyme.missed_cleavages,
                                0..=5,
                            )
                            .text("Missed Cleavages"),
                        );
                        ui.add(
                            egui::Slider::new(&mut self.config.database.enzyme.min_len, 1..=20)
                                .text("Min Length"),
                        );
                        ui.add(
                            egui::Slider::new(&mut self.config.database.enzyme.max_len, 20..=100)
                                .text("Max Length"),
                        );
                        ui.horizontal(|ui| {
                            ui.label("Cleave At:");
                            ui.text_edit_singleline(&mut self.config.database.enzyme.cleave_at);
                        });
                        ui.horizontal(|ui| {
                            ui.label("Restrict:");
                            ui.text_edit_singleline(&mut self.config.database.enzyme.restrict);
                        });
                        ui.checkbox(&mut self.config.database.enzyme.c_terminal, "C-Terminal");
                    });

                    // Mass Ranges
                    ui.group(|ui| {
                        ui.heading("Mass Ranges");
                        ui.add(
                            egui::Slider::new(
                                &mut self.config.database.fragment_min_mz,
                                100.0..=500.0,
                            )
                            .text("Fragment Min m/z"),
                        );
                        ui.add(
                            egui::Slider::new(
                                &mut self.config.database.fragment_max_mz,
                                1000.0..=3000.0,
                            )
                            .text("Fragment Max m/z"),
                        );
                        ui.add(
                            egui::Slider::new(
                                &mut self.config.database.peptide_min_mass,
                                300.0..=1000.0,
                            )
                            .text("Peptide Min Mass"),
                        );
                        ui.add(
                            egui::Slider::new(
                                &mut self.config.database.peptide_max_mass,
                                3000.0..=7000.0,
                            )
                            .text("Peptide Max Mass"),
                        );
                    });

                    ui.checkbox(&mut self.config.database.generate_decoys, "Generate Decoys");
                });

                // Tolerance Configuration Section
                ui.collapsing("Tolerance Settings", |ui| {
                    ui.group(|ui| {
                        ui.heading("Precursor Tolerance");
                        ui.radio_value(
                            &mut self.precursor_tolerance_type,
                            ToleranceType::PPM,
                            "PPM",
                        );
                        ui.radio_value(&mut self.precursor_tolerance_type, ToleranceType::DA, "Da");

                        match self.precursor_tolerance_type {
                            ToleranceType::PPM => {
                                self.config.precursor_tol = ToleranceConfig::Ppm(-10.0, 10.0)
                            }
                            ToleranceType::DA => {
                                self.config.precursor_tol = ToleranceConfig::Da(-0.5, 0.5)
                            }
                        }
                    });

                    ui.group(|ui| {
                        ui.heading("Fragment Tolerance");
                        ui.radio_value(
                            &mut self.fragment_tolerance_type,
                            ToleranceType::PPM,
                            "PPM",
                        );
                        ui.radio_value(&mut self.fragment_tolerance_type, ToleranceType::DA, "Da");

                        match self.fragment_tolerance_type {
                            ToleranceType::PPM => {
                                self.config.fragment_tol = ToleranceConfig::Ppm(-10.0, 10.0)
                            }
                            ToleranceType::DA => {
                                self.config.fragment_tol = ToleranceConfig::Da(-0.5, 0.5)
                            }
                        }
                    });
                });

                // General Settings Section
                ui.collapsing("General Settings", |ui| {
                    ui.add(egui::Slider::new(&mut self.config.min_peaks, 5..=50).text("Min Peaks"));
                    ui.add(
                        egui::Slider::new(&mut self.config.max_peaks, 50..=500).text("Max Peaks"),
                    );
                    ui.add(
                        egui::Slider::new(&mut self.config.min_matched_peaks, 3..=20)
                            .text("Min Matched Peaks"),
                    );
                    ui.add(
                        egui::Slider::new(&mut self.config.max_fragment_charge, 1..=5)
                            .text("Max Fragment Charge"),
                    );
                    ui.add(
                        egui::Slider::new(&mut self.config.report_psms, 1..=10).text("Report PSMs"),
                    );

                    ui.checkbox(&mut self.config.deisotope, "Deisotope");
                    ui.checkbox(&mut self.config.chimera, "Chimera");
                    ui.checkbox(&mut self.config.wide_window, "Wide Window");
                    ui.checkbox(&mut self.config.predict_rt, "Predict RT");
                });

                ui.add_space(20.0);

                ui.horizontal(|ui| {
                    let launch_button = ui.add_enabled(
                        self.process.is_none(), // Disable when process is running
                        egui::Button::new("Launch"),
                    );

                    if launch_button.clicked() {
                        self.status_message = match self.launch_application() {
                            Ok(_) => "Application launched successfully!".to_string(),
                            Err(e) => format!("Error: {}", e),
                        };
                    }
                });

                if !self.status_message.is_empty() {
                    ui.colored_label(
                        if self.status_message.starts_with("Error") {
                            egui::Color32::RED
                        } else {
                            egui::Color32::GREEN
                        },
                        &self.status_message,
                    );
                }
            });
        });

        // Request continuous repaint while process is running
        if self.process.is_some() {
            ctx.request_repaint();
        }
    }
}

impl SageLauncher {
    fn update_process_status(&mut self) {
        if let Some((process, start_time)) = &mut self.process {
            // Check if process is still running
            match process.try_wait() {
                Ok(Some(status)) => {
                    // Process has finished
                    self.status_message = format!("Process finished with status: {}", status);
                    self.process = None;
                    self.elapsed_time.clear();
                }
                Ok(None) => {
                    // Process still running, update elapsed time
                    let elapsed = start_time.elapsed();
                    self.elapsed_time = format_duration(elapsed);
                }
                Err(e) => {
                    // Error checking process status
                    self.status_message = format!("Error checking process status: {}", e);
                    self.process = None;
                    self.elapsed_time.clear();
                }
            }
        }
    }

    fn stop_process(&mut self) {
        if let Some((mut process, _)) = self.process.take() {
            match process.kill() {
                Ok(_) => {
                    self.status_message = "Process stopped".to_string();
                }
                Err(e) => {
                    self.status_message = format!("Error stopping process: {}", e);
                }
            }
            self.elapsed_time.clear();
        }
    }

    fn launch_application(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        // Validate inputs
        if self.sage_executable_path.is_empty() {
            return Err("Sage executable is not selected".into());
        }
        if self.config.database.fasta.is_empty() {
            return Err("FASTA file is not selected".into());
        }
        if self.config.mzml_paths.is_empty() {
            return Err("mzML file is not selected".into());
        }

        // Create config JSON
        let config_json = serde_json::to_string_pretty(&self.config)?;

        // Save config to temporary file
        let config_path = std::env::temp_dir().join("sage_config.json");
        fs::write(&config_path, config_json.clone())?;

        // Print the configuration (for debugging)
        println!("Generated configuration:");
        println!("{}", config_json);

        // Launch the process
        let child = Command::new(&self.sage_executable_path)
            .arg(config_path.to_str().unwrap())
            .spawn()?;

        // Store process handle and start time
        self.process = Some((child, Instant::now()));

        Ok(())
    }
}

// Helper function to format duration
fn format_duration(duration: Duration) -> String {
    let total_secs = duration.as_secs();
    let hours = total_secs / 3600;
    let minutes = (total_secs % 3600) / 60;
    let seconds = total_secs % 60;

    if hours > 0 {
        format!("{}h {}m {}s", hours, minutes, seconds)
    } else if minutes > 0 {
        format!("{}m {}s", minutes, seconds)
    } else {
        format!("{}s", seconds)
    }
}

fn main() -> Result<(), eframe::Error> {
    let options = eframe::NativeOptions {
        // initial_window_size: Some(egui::vec2(600.0, 800.0)),
        ..Default::default()
    };

    eframe::run_native(
        "Sage Launcher",
        options,
        Box::new(|_cc| Ok(Box::new(SageLauncher::default()))),
    )
}
