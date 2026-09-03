"""
multiTargetGUI.py

Desktop GUI (Tkinter) for multiTargetLAFCA.py.

Lets you set:
    - number of targets
    - each target's ratio (comma-separated fractions)
    - error tolerance (err)
    - animation fps

then runs the real pipeline (Z3 tree generation -> LAFCA placement ->
parallel DFL routing -> animation) in a BACKGROUND THREAD so the window
stays responsive, streams the pipeline's print() output into a log box,
and plays the resulting sim.gif directly in the window once done.

Requirements: multiTargetLAFCA.py, animateMultiTarget.py, and everything
they depend on (skewedTreeGenerator, createTreeForShared, getLoadingData,
LAFCADFL package, z3-solver) must be on the Python path -- run this file
from the same folder (mtp_july) as those modules.
Pillow (PIL) is required for GIF playback -- it's already a dependency
of matplotlib's PillowWriter, used to SAVE the gif, so it should already
be installed. If not: pip install Pillow
"""

import os
import sys
import time
import queue
import threading
import tempfile
import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox

from PIL import Image, ImageTk

import multiTargetLAFCA as mtl


class _RedirectText:
    """File-like object that pushes print() output into a thread-safe queue."""
    def __init__(self, q):
        self.q = q

    def write(self, s):
        if s.strip():
            self.q.put(("log", s))

    def flush(self):
        pass


class MultiTargetGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Multi-Target Skewed Tree Simulator")
        self.geometry("1100x900")
        self.minsize(800, 650)

        self.queue = queue.Queue()
        self.ratio_entries = []
        self.gif_pil_frames = []
        self.gif_durations = []
        self.gif_index = 0
        self.gif_job = None

        self._build_controls()
        self._build_log()
        self._build_gif_view()

        self.after(150, self._poll_queue)

    # ----- UI construction ---------------------------------------------
    def _build_controls(self):
        frm = ttk.LabelFrame(self, text="Setup")
        frm.pack(side="top", fill="x", padx=10, pady=8)

        ttk.Label(frm, text="Number of targets:").grid(row=0, column=0, sticky="w", padx=4, pady=4)
        self.k_var = tk.IntVar(value=3)
        ttk.Spinbox(frm, from_=1, to=9, width=5, textvariable=self.k_var,
                    command=self._rebuild_ratio_rows).grid(row=0, column=1, sticky="w", padx=4, pady=4)

        ttk.Label(frm, text="Error (err):").grid(row=0, column=2, sticky="w", padx=4, pady=4)
        self.err_var = tk.StringVar(value="0.01")
        ttk.Entry(frm, textvariable=self.err_var, width=8).grid(row=0, column=3, sticky="w", padx=4, pady=4)

        ttk.Label(frm, text="Animation fps:").grid(row=0, column=4, sticky="w", padx=4, pady=4)
        self.fps_var = tk.StringVar(value="2")
        ttk.Entry(frm, textvariable=self.fps_var, width=5).grid(row=0, column=5, sticky="w", padx=4, pady=4)

        ttk.Label(frm, text="Playback speed:").grid(row=1, column=0, sticky="w", padx=4, pady=4)
        self.speed_var = tk.DoubleVar(value=1.0)
        speed_scale = ttk.Scale(frm, from_=0.1, to=3.0, orient="horizontal",
                                 variable=self.speed_var, length=160,
                                 command=lambda v: self.speed_label.configure(
                                     text=f"{float(v):.2f}x"))
        speed_scale.grid(row=1, column=1, columnspan=2, sticky="w", padx=4, pady=4)
        self.speed_label = ttk.Label(frm, text="1.00x")
        self.speed_label.grid(row=1, column=3, sticky="w", padx=4, pady=4)

        self.ratio_frame = ttk.LabelFrame(
            self, text="Target ratios (comma-separated fractions, e.g. 0.5,0.25,0.25)")
        self.ratio_frame.pack(side="top", fill="x", padx=10, pady=4)
        self._rebuild_ratio_rows()

        btn_frame = ttk.Frame(self)
        btn_frame.pack(side="top", fill="x", padx=10, pady=6)
        self.run_btn = ttk.Button(btn_frame, text="Run Simulation", command=self._on_run)
        self.run_btn.pack(side="left")
        self.status_var = tk.StringVar(value="Idle")
        ttk.Label(btn_frame, textvariable=self.status_var).pack(side="left", padx=12)

    def _rebuild_ratio_rows(self):
        for child in self.ratio_frame.winfo_children():
            child.destroy()
        self.ratio_entries = []
        k = self.k_var.get()
        defaults = ["0.5,0.25,0.25", "0.25,0.5,0.25", "0.125,0.375,0.5"]
        for i in range(k):
            ttk.Label(self.ratio_frame, text=f"Target {i + 1}:").grid(
                row=i, column=0, sticky="w", padx=4, pady=2)
            var = tk.StringVar(value=defaults[i] if i < len(defaults) else "")
            ttk.Entry(self.ratio_frame, textvariable=var, width=45).grid(
                row=i, column=1, sticky="w", padx=4, pady=2)
            self.ratio_entries.append(var)

    def _build_log(self):
        frm = ttk.LabelFrame(self, text="Log")
        frm.pack(side="top", fill="both", expand=False, padx=10, pady=4)
        self.log_box = scrolledtext.ScrolledText(frm, height=7, state="disabled",
                                                   font=("Consolas", 9))
        self.log_box.pack(fill="both", expand=True, padx=4, pady=4)

    def _build_gif_view(self):
        self.warning_var = tk.StringVar(value="")
        self.warning_label = tk.Label(self, textvariable=self.warning_var,
                                       bg="#fff3cd", fg="#7a5b00",
                                       font=("Segoe UI", 9, "bold"),
                                       anchor="w", padx=8, pady=4)
        # packed/unpacked dynamically -- only shown when there's something to warn about

        notebook = ttk.Notebook(self)
        notebook.pack(side="top", fill="both", expand=True, padx=10, pady=4)
        self.notebook = notebook

        sim_frame = ttk.Frame(notebook)
        notebook.add(sim_frame, text="Simulation")
        self.gif_label = ttk.Label(sim_frame, text="(simulation will appear here after a run)",
                                    anchor="center")
        self.gif_label.pack(fill="both", expand=True, padx=4, pady=4)

        trees_frame = ttk.Frame(notebook)
        notebook.add(trees_frame, text="Trees")
        self.trees_canvas_frame = trees_frame
        self.tree_photo_refs = []   # keep PhotoImage references alive
        ttk.Label(trees_frame, text="(individual target tree diagrams will appear here after a run)",
                  anchor="center").pack(fill="both", expand=True, padx=4, pady=4)

    def _show_warning(self, text):
        self.warning_var.set(text)
        self.warning_label.pack(side="top", fill="x", before=self.notebook)

    def _clear_warning(self):
        self.warning_var.set("")
        self.warning_label.pack_forget()

    def _show_trees(self, tree_image_paths):
        for child in self.trees_canvas_frame.winfo_children():
            child.destroy()
        self.tree_photo_refs = []

        if not tree_image_paths:
            ttk.Label(self.trees_canvas_frame, text="(no tree diagrams available)",
                      anchor="center").pack(fill="both", expand=True, padx=4, pady=4)
            return

        # scrollable strip of tree thumbnails -- BOTH directions, since a
        # deep tree's rendered PNG (graphviz, rankdir='BT') can easily be
        # taller than the available panel height, and previously there was
        # no vertical scrollbar at all, so anything below the visible area
        # was simply unreachable, not just visually cropped.
        outer = ttk.Frame(self.trees_canvas_frame)
        outer.pack(side="top", fill="both", expand=True)
        outer.rowconfigure(0, weight=1)
        outer.columnconfigure(0, weight=1)

        canvas = tk.Canvas(outer, highlightthickness=0)
        vbar = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        hbar = ttk.Scrollbar(outer, orient="horizontal", command=canvas.xview)
        canvas.configure(yscrollcommand=vbar.set, xscrollcommand=hbar.set)

        canvas.grid(row=0, column=0, sticky="nsew")
        vbar.grid(row=0, column=1, sticky="ns")
        hbar.grid(row=1, column=0, sticky="ew")

        # mouse wheel support, scoped to only when hovering this canvas --
        # bind_all would leak a new global binding (referencing an
        # already-destroyed canvas) every time a new run replaces this tab
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        def _bind_wheel(_e):
            canvas.bind_all("<MouseWheel>", _on_mousewheel)
        def _unbind_wheel(_e):
            canvas.unbind_all("<MouseWheel>")
        canvas.bind("<Enter>", _bind_wheel)
        canvas.bind("<Leave>", _unbind_wheel)

        inner = ttk.Frame(canvas)
        canvas.create_window((0, 0), window=inner, anchor="nw")

        # conservative default cap -- still scrollable if a tree is taller
        # than this, but keeps the common case fully visible without needing
        # to scroll at all
        THUMB_MAX = (300, 420)

        for idx, path in enumerate(tree_image_paths):
            col = ttk.Frame(inner, padding=6)
            col.pack(side="left", fill="y", anchor="n")
            ttk.Label(col, text=f"tr{idx + 1}", font=("Segoe UI", 9, "bold")).pack(side="top")
            try:
                img = Image.open(path)
                img.thumbnail(THUMB_MAX, Image.LANCZOS)
                photo = ImageTk.PhotoImage(img)
                self.tree_photo_refs.append(photo)
                ttk.Label(col, image=photo).pack(side="top")
            except Exception as e:
                ttk.Label(col, text=f"(couldn't load: {e})").pack(side="top")

        inner.update_idletasks()
        canvas.configure(scrollregion=canvas.bbox("all"))

    # ----- logging -------------------------------------------------------
    def _log(self, msg):
        self.log_box.configure(state="normal")
        self.log_box.insert("end", msg if msg.endswith("\n") else msg + "\n")
        self.log_box.see("end")
        self.log_box.configure(state="disabled")

    # ----- run -------------------------------------------------------------
    def _on_run(self):
        ratios = []
        for i, var in enumerate(self.ratio_entries):
            text = var.get().strip()
            if not text:
                messagebox.showerror("Missing input", f"Target {i + 1}'s ratio is empty.")
                return
            try:
                vals = [float(v) for v in text.split(",") if v.strip() != ""]
            except ValueError:
                messagebox.showerror("Invalid input",
                                      f"Target {i + 1}'s ratio isn't valid comma-separated numbers.")
                return
            if len(vals) < 2:
                messagebox.showerror("Invalid input",
                                      f"Target {i + 1} needs at least 2 reagent fractions.")
                return
            ratios.append(vals)

        try:
            err = float(self.err_var.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Error must be a number, e.g. 0.01.")
            return

        try:
            fps = int(self.fps_var.get())
        except ValueError:
            fps = 2

        self.run_btn.configure(state="disabled")
        self.status_var.set("Running...")
        self._stop_gif()
        self._clear_warning()
        self.gif_label.configure(image="", text="Running simulation...")
        self.log_box.configure(state="normal")
        self.log_box.delete("1.0", "end")
        self.log_box.configure(state="disabled")

        thread = threading.Thread(target=self._run_pipeline, args=(ratios, err, fps), daemon=True)
        thread.start()

    def _run_pipeline(self, ratios, err, fps):
        old_stdout = sys.stdout
        sys.stdout = _RedirectText(self.queue)
        try:
            tmpdir = tempfile.mkdtemp(prefix="mtl_gui_")
            inputFiles = []
            for i, vals in enumerate(ratios):
                path = os.path.join(tmpdir, f"input{i + 1}.csv")
                with open(path, "w") as f:
                    f.write(",".join(str(v) for v in vals) + "\n")
                inputFiles.append(path)

            batch = f"run_{int(time.time())}"
            directory = "./gui_output/"
            result = mtl.multiTargetLAFCA(
                inputFiles=tuple(inputFiles), err=err,
                directory=directory, batch=batch,
                animate=True, fps=fps,
            )

            if result is None:
                self.queue.put(("error", "No SAT tree could be generated for any target -- "
                                          "check the log above for details."))
            else:
                self.queue.put(("done", result))
        except Exception as e:
            self.queue.put(("error", f"{type(e).__name__}: {e}"))
        finally:
            sys.stdout = old_stdout

    # ----- queue polling (main thread only touches the UI) -----------------
    def _poll_queue(self):
        try:
            while True:
                item = self.queue.get_nowait()
                kind = item[0]
                if kind == "log":
                    self._log(item[1])
                elif kind == "done":
                    result = item[1]
                    self.status_var.set(
                        f"Done -- Flow={result['K']} Bendings={result['B']} PathLength={result['L']}")
                    self.run_btn.configure(state="normal")

                    self._show_trees(result.get("treeImages", []))

                    skipped = result.get("skipped", [])
                    if skipped:
                        names = ", ".join(s["path"] for s in skipped)
                        self._show_warning(
                            f"\u26a0 {len(skipped)} of {len(self.ratio_entries)} target(s) could NOT "
                            f"be synthesized (no valid tree found) and were SKIPPED: {names}")
                    else:
                        self._clear_warning()

                    if result.get("simPath"):
                        self._load_gif(result["simPath"])
                elif kind == "error":
                    self.status_var.set("Error")
                    self.run_btn.configure(state="normal")
                    self._log(f"ERROR: {item[1]}")
                    messagebox.showerror("Simulation failed", item[1])
        except queue.Empty:
            pass
        self.after(150, self._poll_queue)

    # ----- gif playback ------------------------------------------------------
    def _load_gif(self, path):
        if not os.path.exists(path):
            self._log(f"animation file not found: {path}")
            self.gif_label.configure(text=f"(animation file not found: {path})")
            return
        img = Image.open(path)
        pil_frames = []
        durations = []
        try:
            while True:
                pil_frames.append(img.copy().convert("RGBA"))
                durations.append(img.info.get("duration", 200))  # ms, as saved
                img.seek(img.tell() + 1)
        except EOFError:
            pass
        self.gif_pil_frames = pil_frames
        self.gif_durations = durations
        self.gif_index = 0
        self.gif_label.configure(text="")
        self._animate_gif()

    def _fit_image(self, pil_img):
        """Scale (up or down) to fill the label's CURRENT allocated space,
        preserving aspect ratio, so nothing is cropped by the widget
        boundary regardless of window size."""
        w = self.gif_label.winfo_width()
        h = self.gif_label.winfo_height()
        if w < 10 or h < 10:   # not laid out yet on the very first frame
            w, h = 900, 520
        iw, ih = pil_img.size
        scale = min(w / iw, h / ih)
        new_size = (max(1, int(iw * scale)), max(1, int(ih * scale)))
        return pil_img.resize(new_size, Image.LANCZOS)

    def _animate_gif(self):
        if not self.gif_pil_frames:
            return
        fitted = self._fit_image(self.gif_pil_frames[self.gif_index])
        photo = ImageTk.PhotoImage(fitted)
        self.gif_label.configure(image=photo)
        self.gif_label.image = photo   # keep a reference, tkinter needs this

        base_delay = self.gif_durations[self.gif_index]
        speed = max(self.speed_var.get(), 0.05)   # guard against 0/negative
        delay = max(20, int(base_delay / speed))  # 20ms floor, avoids runaway loop

        self.gif_index = (self.gif_index + 1) % len(self.gif_pil_frames)
        self.gif_job = self.after(delay, self._animate_gif)

    def _stop_gif(self):
        if self.gif_job is not None:
            self.after_cancel(self.gif_job)
            self.gif_job = None
        self.gif_pil_frames = []
        self.gif_durations = []


if __name__ == "__main__":
    app = MultiTargetGUI()
    app.mainloop()