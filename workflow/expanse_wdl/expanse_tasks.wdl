version 1.0

task extract_field {
  input {
    String script_dir
    File script = "~{script_dir}/main_dataset/decompress_trait.py"
    File ukbconv = "~{script_dir}/ukb_utilities/ukbconv"
    File encoding = "~{script_dir}/ukb_utilities/encoding.ukb"
    Array[File]+ fields_files = ["~{script_dir}/main_dataset/raw_data/fields46781.ukb", "~{script_dir}/main_dataset/raw_data/fields46782.ukb"]
    Array[File]+ enc_files = ["~{script_dir}/main_dataset/raw_data/ukb46781.enc_ukb", "~{script_dir}/main_dataset/raw_data/ukb46782.enc_ukb"]

    Int id # data field id
  }

  output {
    File data = "~{id}_cleaned.tab"
  }

  command <<<
    envsetup ~{script} \
      ~{id} \
      ~{id} \
      ~{ukbconv} \
      ~{encoding} \
      --fields-files ~{sep=" " fields_files} \
      --enc-files ~{sep=" " enc_files} && \
    sed -e 's/^\t//' -e 's/\t$//' ~{id}.txt > ~{id}_cleaned.tab
    # ukb's extract skip pad every line but the first with an extra tab at the beginning and end
    # seemingly for no reason. This messes up tsv column alignment, so remove them.
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "5h"
    memory: "2GB"
  }
}

task rename_file {
  input {
    File f
    String name
  }

  output {
    File out = name
  }

  command <<<
    ln ~{f} ~{name}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "1h"
    memory: "2GB"
    shortTask: true
  }
}
