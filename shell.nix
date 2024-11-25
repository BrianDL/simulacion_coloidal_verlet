{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = with pkgs; [
    entr
    zig
    python312
    python312Packages.pip
    stdenv.cc.cc.lib  # This includes libstdc++
    zlib  # Add zlib to resolve libz.so.1 issue
  ];

  shellHook = ''
    # Create a virtual environment if it doesn't exist
    if [ ! -d "venv" ]; then
      python3 -m venv .venv
    fi

    # Activate the virtual environment
    source .venv/bin/activate

    # Upgrade pip
    pip install --upgrade pip

    # Install requirements if requirements.txt exists
    if [ -f requirements.txt ]; then
      pip install -r requirements.txt
    else
      echo "requirements.txt not found. Please create one with your project dependencies."
    fi

    # Set LD_LIBRARY_PATH to include libstdc++ and zlib
    export LD_LIBRARY_PATH=${pkgs.stdenv.cc.cc.lib}/lib:${pkgs.zlib}/lib:$LD_LIBRARY_PATH

    # Inform the user that the environment is ready
    echo "Python virtual environment is set up and requirements are installed."
    echo "To deactivate the virtual environment, type 'deactivate'."
  '';
}